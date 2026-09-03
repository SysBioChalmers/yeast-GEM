function model = commitYeastModel(model,formats,allowNoGrowth)
% commitYeastModel
%   Independently re-applies the minimal_Y6 canonical medium, SBO terms
%   and aerobic/anaerobic growth checks (does not trust that
%   saveYeastYaml already ran), validates the model as SBML, and writes
%   the requested file formats. Use this for local binary builds or a
%   release; for routine curation saves, use saveYeastYaml instead, which
%   is faster and does not require SBML validation.
%
%   Does not write model/yeast-GEM.yml or the annotation tsvs -- that is
%   exclusively saveYeastYaml's job.
%
% Inputs:
%   model           (struct) model to commit. Preferably RAVEN format,
%                   although COBRA format is also allowed, but some
%                   fields might be lost in the conversion.
%   formats         (cell array, opt) which file formats to write, any of
%                   'xml', 'txt', 'xlsx', 'mat' (default {'xml','txt'}).
%   allowNoGrowth   (bool, opt) if committing should be allowed whenever
%                   the model cannot grow, returning a warning (default
%                   true), otherwise will error.
%
% Output:
%   model   the model, with minimal_Y6 and SBO terms applied.
%
% Usage: model = commitYeastModel(model,formats,allowNoGrowth)

if nargin < 2 || isempty(formats)
    formats = {'xml','txt'};
end
if nargin < 3 || isempty(allowNoGrowth)
    allowNoGrowth = true;
end
if ~(exist('ravenCobraWrapper.m','file')==2)
    error(['RAVEN cannot be found. See README.md for installation '...
        'instructions. RAVEN is required to make sure that the model '...
        'is stored in the correct file formats for use in the '...
        'yeast-GEM GitHub repository'])
end

% Export as RAVEN format
if isfield(model,'rules')
    model = ravenCobraWrapper(model);
end

%Get and change to the script folder, as all folders are relative to this
%folder
scriptFolder = fileparts(which(mfilename));
currentDir = cd(scriptFolder);

%Set minimal media
cd modelCuration
model = minimal_Y6(model);
cd ..

%Update SBO terms in model:
cd missingFields
model = addSBOterms(model);
cd ..

%Check if model is a valid SBML structure:
exportModel(model,'tempModel.xml',false,false,true);
try
    [~,~,errors] = evalc('TranslateSBML_RAVEN(''tempModel.xml'',1,0)');
catch
    [~,~,errors] = evalc('TranslateSBML(''tempModel.xml'',1,0)');
end
if any(strcmp({errors.severity},'Error'))
    delete('tempModel.xml');
    error('Model should be a valid SBML structure. Please fix all errors before committing.')
end

%Check if model can grow:
checkGrowth(model,'aerobic',allowNoGrowth)
checkGrowth(model,'anaerobic',allowNoGrowth)

%Write the requested file formats. .xml is written from the
%already-validated temp file; the others are written directly from model.
if ismember('xml',formats)
    copyfile('tempModel.xml','../model/yeast-GEM.xml')
end
delete('tempModel.xml');

otherFormats = setdiff(formats,{'xml'});
if ~isempty(otherFormats)
    exportForGit(model,'yeast-GEM','../model',otherFormats,false,false);
end

%Write deltaG fields to file
cd missingFields
saveDeltaG(model,false);
cd ..

%Convert notation "e-005" to "e-05 " in stoich. coeffs. to avoid
%inconsistencies between Windows and MAC (.xml only):
if ismember('xml',formats)
    copyfile('../model/yeast-GEM.xml','backup.xml')
    fin  = fopen('backup.xml','r');
    fout = fopen('../model/yeast-GEM.xml','w');
    still_reading = true;
    while still_reading
        inline = fgets(fin);
        if ~ischar(inline)
            still_reading = false;
        else
            if ~isempty(regexp(inline,'[0-9]e-?00[0-9]','once'))
                inline = regexprep(inline,'(?<=[0-9]e-?)00(?=[0-9])','0');
            end
            fwrite(fout,inline);
        end
    end
    fclose('all');
    delete('backup.xml');
end

%Switch back to original folder
cd(currentDir)
end

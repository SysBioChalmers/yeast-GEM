function saveYeastModel(model,upDATE,allowNoGrowth,binaryFiles)
% saveYeastModel
%   Saves model as a .xml, .txt and .yml file. Also updates complementary
%   files (boundaryMets.txt, README.md and dependencies.txt).
%
%   Inputs: model           (struct) model to save. Preferably RAVEN
%                           format, although COBRA format is also allowed,
%                           but some fields might be lost in the
%                           conversion.
%           upDATE          (bool, opt) If updating the date in the README file
%                           is needed (default true)
%           allowNoGrowth   (bool, opt) if saving should be allowed whenever
%                           the model cannot grow, returning a warning (default
%                           = true), otherwise will error
%           binaryFiles     (bool, opt) if the model should be stored in
%                           binary file formats (= xlsx and mat)
%
%   Usage: saveYeastModel(model,upDATE,allowNoGrowth,binaryFiles)
%

if nargin < 2
    upDATE = true;
end
if nargin < 3
    allowNoGrowth = true;
end
if nargin < 4
    binaryFiles = false;
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
[~,~,errors] = evalc('TranslateSBML(''tempModel.xml'',1,0)');
if any(strcmp({errors.severity},'Error'))
    delete('tempModel.xml');
    error('Model should be a valid SBML structure. Please fix all errors before saving.')
end

%Check if model can grow:
checkGrowth(model,'aerobic',allowNoGrowth)
checkGrowth(model,'anaerobic',allowNoGrowth)

%Update .xml, .txt and .yml models:
copyfile('tempModel.xml','../model/yeast-GEM.xml')
delete('tempModel.xml');
if binaryFiles==false
    exportForGit(model,'yeast-GEM','../model',{'yml','txt'},false,false);
else
    exportForGit(model,'yeast-GEM','../model',{'yml','txt','xlsx','mat'},false,false);
end

%Write deltaG fields to file
cd missingFields
saveDeltaG(model,false);
cd ..

%Update README file: date + size of model
modelVersion = regexprep(model.id,'yeastGEM_v?','');
nGenes=num2str(numel(model.genes));
nMets=num2str(numel(model.mets));
nRxns=num2str(numel(model.rxns));
copyfile('../README.md','backup.md')
fin  = fopen('backup.md','r');
fout = fopen('../README.md','w');
newStats = ['| $1 | ' datestr(now,'dd-mmm-yyyy') ' | ' modelVersion ' | ' nRxns ' | ' nMets ' | ' nGenes ' |'];
searchStats = '^\| (\_Saccharomyces cerevisiae\_) \| \d{2}-\D+-\d{4} \| (\d+\.\d+\.\d+|develop) \| \d+ \| \d+ \| \d+ \|';
while ~feof(fin)
    str = fgets(fin);
    inline = regexprep(str,searchStats,newStats);
    inline = unicode2native(inline,'UTF-8');
    fwrite(fout,inline);
end
fclose('all');
delete('backup.md');

%Convert notation "e-005" to "e-05 " in stoich. coeffs. to avoid
%inconsistencies between Windows and MAC:
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

%Switch back to original folder
cd(currentDir)
end

%%
function checkGrowth(model,condition,allowNoGrowth)
%Function that checks if the model can grow or not using COBRA under a
%given condition (aerobic or anaerobic). Will either return warnings or
%errors depending on allowNoGrowth.

if strcmp(condition,'anaerobic')
    cd otherChanges
    model = anaerobicModel(model);
    cd ..
end
try
    xPos = strcmp(model.rxnNames,'growth');
    sol  = solveLP(model);
    if sol.x(xPos) < 1e-6
        dispText = ['The model is not able to support growth under ' ...
            condition ' conditions. Please ensure the model can grow'];
    end
catch
    dispText = ['The model yields an infeasible simulation using COBRA ' ...
        'under ' condition ' conditions. Please ensure the model ' ...
        'can be simulated with COBRA'];
end

if exist('dispText','var')
    if allowNoGrowth
        warning([dispText ' before opening a PR.'])
    else
        error([dispText ' before committing.'])
    end
end
end

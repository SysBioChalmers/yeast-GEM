%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saveYeastModel(model,upDATE)
% Saves model as a .xml, .txt and .yml file. Also updates complementary
% files (boundaryMets.txt, README.md and dependencies.txt).
%
% model     model structure to save (note: must be in COBRA format)
% upDATE    logical =true if updating the date in the README file is needed
%           (opt, default true)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function saveYeastModel(model,upDATE)

if nargin < 2
    upDATE = true;
end

%Get and change to the script folder, as all folders are relative to this
%folder
scriptFolder = fileparts(which(mfilename));
currentDir = cd(scriptFolder);

%Set minimal media
cd modelCuration
model = minimal_Y6(model);
cd ..

%Delete model.grRules (redundant and possibly conflicting with model.rules):
if isfield(model,'grRules')
    model = rmfield(model,'grRules');
end

%Update SBO terms in model:
cd missingFields
model = addSBOterms(model);
cd ..

%Save "proteins" ("fbc:name" in the xml file) = "geneNames" ("fbc:label" in the xml file):
model.proteins = model.geneNames;

%Check if model is a valid SBML structure:
writeCbModel(model,'sbml','tempModel.xml');
[~,errors] = TranslateSBML('tempModel.xml');
if ~isempty(errors)
    delete('tempModel.xml');
    error('Model should be a valid SBML structure. Please fix all errors before saving.')
end

%Update .xml, .txt and .yml models:
copyfile('tempModel.xml','../ModelFiles/xml/yeastGEM.xml')
delete('tempModel.xml');
writeCbModel(model,'text','../ModelFiles/txt/yeastGEM.txt');
exportForGit(model,'yeastGEM','..',{'yml'});

%Detect boundary metabolites and save them in a .txt file:
fid = fopen('../ModelFiles/boundaryMets.txt','wt');
for i = 1:length(model.rxns)
    pos = find(model.S(:,i) ~= 0);
    if length(pos) == 1 %Exchange rxn
        fprintf(fid,[model.mets{pos} '\t' model.metNames{pos} '\n']);
    end
end
fclose(fid);

%Update README file: date + size of model
copyfile('../README.md','backup.md')
fin  = fopen('backup.md','r');
fout = fopen('../README.md','w');
still_reading = true;
while still_reading
    inline = fgets(fin);
    if ~ischar(inline)
        still_reading = false;
    else
        if startsWith(inline,'* Last update: ') && upDATE
            inline = ['* Last update: ' datestr(datetime,'yyyy-mm-dd') newline];
        elseif startsWith(inline,'|_Saccharomyces cerevisiae_|')
            inline = ['|_Saccharomyces cerevisiae_|[Yeast 7.6]' ...
                '(https://sourceforge.net/projects/yeast/)|' ...
                num2str(length(model.rxns)) '|' ...
                num2str(length(model.mets)) '|' ...
                num2str(length(model.genes)) '|' newline];
        end
        fwrite(fout,inline);
    end
end
fclose('all');
delete('backup.md');

%Convert notation "e-005" to "e-05 " in stoich. coeffs. to avoid
%inconsistencies between Windows and MAC:
copyfile('../ModelFiles/xml/yeastGEM.xml','backup.xml')
fin  = fopen('backup.xml','r');
fout = fopen('../ModelFiles/xml/yeastGEM.xml','w');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

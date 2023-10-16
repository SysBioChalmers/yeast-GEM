% This scripts applies curations to be applied on yeast-GEM release 8.7.0, to
% get to yeast-GEM release 8.7.1.
% Indicate which Issue/PR are addressed. If multiple curations are performed
% before a new release is made, just add the required code to this script. If
% more extensive coding is required, you can write a separate (generic) function
% that can be kept in the /code/modelCuration folder. Otherwise, try to use
% existing functions whenever possible. In particular /code/curateMetsRxnsGenes
% can do many types of curation.

%% Load yeast-GEM 8.7.0 (requires local yeast-GEM git repository)
cd ..
codeDir=pwd();
model = getEarlierModelVersion('8.7.0');
model.id='yeastGEM_develop';
dataDir=fullfile(pwd(),'..','data','modelCuration','v8_7_1');
cd modelCuration

%% Add Uniprot IDs (PR #349)
% Gather Uniprot IDs with get_uniprot_id.py that is in dataDir.
uniprotFile = fullfile(dataDir,'SGD_with_Uniprot.csv');
fid = fopen(uniprotFile);
uniprot = textscan(fid,'%q %q %q','Delimiter',',','HeaderLines',1);
fclose(fid);

% Only keep data for genes that are in the model
uniprotGenes = uniprot{1,1};
uniprotIDs   = uniprot{1,3};
isInModel    = ismember(uniprotGenes,model.genes);
model        = editMiriam(model,'gene',uniprotGenes(isInModel),'uniprot',uniprotIDs(isInModel),'add');

%% Corrent reaction name r_1024
rxnIdx = getIndexes(model,'r_1024','rxns');
model.rxnNames{rxnIdx} = 'sucrose hydrolyzing enzyme';

%% DO NOT CHANGE OR REMOVE THE CODE BELOW THIS LINE.
% Show some metrics:
cd(fullfile(codeDir,'modelTests'))
disp('Run gene essentiality analysis')
[new.accuracy,new.tp,new.tn,new.fn,new.fp] = essentialGenes(model);
fprintf('Genes in model: %d\n',numel(model.genes));
fprintf('Gene essentiality accuracy: %.4f\n', new.accuracy);
fprintf('True non-essential genes: %d\n', numel(new.tp));
fprintf('True essential genes: %d\n', numel(new.tn));
fprintf('False non-essential genes: %d\n', numel(new.fp));
fprintf('False essential genes: %d\n', numel(new.fn));
fprintf('\nRun growth analysis\n')
R2=growth(model);
fprintf('R2 of growth prediction: %.4f\n', R2);

% Save model:
cd ..
saveYeastModel(model)
cd modelCuration

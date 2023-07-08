% This scripts applies curations to be applied on yeast-GEM release develop.
% Indicate which Issue/PR are addressed. If multiple curations are performed
% before a new release is made, just add the required code to this script. If
% more extensive coding is required, you can write a separate (generic) function
% that can be kept in the /code/modelCuration folder. Otherwise, try to use
% existing functions whenever possible. In particular /code/curateMetsRxnsGenes
% can do many types of curation.

%% Load yeast-GEM develop (requires local yeast-GEM git repository)
codeDir = pwd;
cd ..
model = getEarlierModelVersion('662cdefcb2e8af2594f8c30ce62a75ef6f4ef0de');
model.id='yeastGEM_develop';
dataDir=fullfile(pwd(),'..','data','modelCuration','v8_6_4');
cd modelCuration

%% Add deltaG for reactions and metabolites (PR #330)
cd(dataDir)
%% Gather yETFL data
% First-time run:
websave('input_models.zip','https://zenodo.org/record/4778047/files/input_models.zip?download=1');
unzip('input_models.zip');
modelYETFL = load('input_models/yeast8_thermo_curated.mat');
modelYETFL = modelYETFL.model;
modelYETFL.metNames = regexprep(modelYETFL.metNames,' \[[\w ]+\]$','');
[modelYETFL.metNames,idx]=unique(modelYETFL.metNames(1:end-1));
modelYETFL.metDeltaGFtr = modelYETFL.metDeltaGFtr(idx);
yetfl_metG = array2table([modelYETFL.metNames, num2cell(modelYETFL.metDeltaGFtr)]);
writetable(yetfl_metG,'yetfl_metG.csv');
yetfl_rxnG = array2table([modelYETFL.rxns, num2cell(modelYETFL.rxnDeltaGR)]);
writetable(yetfl_rxnG,'yetfl_rxnG.csv');
clear modelYETFL idx   
rmdir input_models s
delete input_models.zip
% Repeated runs can reuse this file, but will not submit to GitHub
yetfl_metG = readtable('yetfl_metG.csv');
yetfl_rxnG = readtable('yetfl_rxnG.csv');

%% Gather ModelSEED data via get_seed_data.py, CSV obtained how?
seed_metG = readtable('modelseed_metG.csv');

%% Gather dG-predictor data via run_dgpredictor.py
dgpred_rxnG = readtable('dgpred_rxnG.csv');

%% Assign metDeltaG, preferred source: yETFL > ModelSEED
model.metDeltaG = nan(numel(model.mets),1);

% Map to yETFL by metNames
[a,b] = ismember(model.metNames,yetfl_metG.Var1);
model.metDeltaG(a) = round(yetfl_metG.Var2(b(a)),2);

% If no deltaG is assigned, map to ModelSEED by KEGG IDs
noDeltaG = isnan(model.metDeltaG);

kegg = extractMiriam(model.metMiriams,'kegg.compound');
hasKEGG = ~cellfun(@isempty,kegg);
toCheck = find(hasKEGG & noDeltaG);

[a,b] = ismember(kegg(toCheck),seed_metG.kegg);
model.metDeltaG(toCheck(a)) = round(seed_metG.deltag(b(a)),2);

%% Assign rxnDeltaG, preferred source: yETFL > dG-predictor
model.rxnDeltaG = nan(numel(model.rxns),1);

% Map to yETFL by reaction IDs
[a,b] = ismember(model.rxns,yetfl_rxnG.Var1);
model.rxnDeltaG(a) = round(yetfl_rxnG.Var2(b(a)),2);

% If no deltaG is assigned, calculate by dG-predictor
noDeltaG = isnan(model.rxnDeltaG);
% Actually, no deltaG is found for model.rxns(3983:4131). Current
% dG-predictor dataset has data for model.rxns(1:4062). For now, use the
% deltaGs for model.rxns(3983:4062).
model.rxnDeltaG(3983:4062) = round(dgpred_rxnG.detaG(3983:4062),2);

%% RAVEN does not yet support .rxnDeltaG and .metDeltaG fields.
% For now, write these two tables. RAVEN will soon support I/O of deltaG in
% YAML files, and a custom function could be written to add the field after
% loading SBML file as well. Otherwise, the current script can be run
cd ../../databases/
metG = array2table([model.mets, num2cell(model.metDeltaG)]);
writetable(metG,'model_metDeltaG.csv');
rxnG = array2table([model.rxns, num2cell(model.rxnDeltaG)]);
writetable(rxnG,'model_rxnDeltaG.csv');
cd(dataDir)

% %% DO NOT CHANGE OR REMOVE THE CODE BELOW THIS LINE.
% % Show some metrics:
% cd ../modelTests
% disp('Run gene essentiality analysis')
% [new.accuracy,new.tp,new.tn,new.fn,new.fp] = essentialGenes(model);
% fprintf('Genes in model: %d\n',numel(model.genes));
% fprintf('Gene essentiality accuracy: %.4f\n', new.accuracy);
% fprintf('True non-essential genes: %d\n', numel(new.tp));
% fprintf('True essential genes: %d\n', numel(new.tn));
% fprintf('False non-essential genes: %d\n', numel(new.fp));
% fprintf('False essential genes: %d\n', numel(new.fn));
% fprintf('\nRun growth analysis\n')
% R2=growth(model);
% fprintf('R2 of growth prediction: %.4f\n', R2);
% 
% % Save model:
cd(fullfile(codeDir,'..'))
saveYeastModel(model)
cd modelCuration

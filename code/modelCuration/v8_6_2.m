% This scripts applies curations to be applied on yeast-GEM release 8.6.2.
% Indicate which Issue/PR are addressed. If multiple curations are performed
% before a new release is made, just add the required code to this script. If
% more extensive coding is required, you can write a separate (generic) function
% that can be kept in the /code/modelCuration folder. Otherwise, try to use
% existing functions whenever possible. In particular /code/curateMetsRxnsGenes
% can do many types of curation.

%% Load yeast-GEM 8.6.2 (requires local yeast-GEM git repository)
cd ..
model = getEarlierModelVersion('8.6.2');
model.id='yeastGEM_develop';
dataDir=fullfile(pwd(),'..','data','modelCuration','v8_6_2');
cd modelCuration

%% Volatile Esters & Polyphosphate Reactions (PR #336)
% See https://github.com/SysBioChalmers/yeast-GEM/pull/336 for more detailed
% explanation of what changes were made by including 8 new distinct ester
% reactions and  polyphosphate synthesis and transport reactions.

%% Curate gene association for transport rxns (PR #306)
% Add new reactions and genes
metsInfo = fullfile(dataDir,'VolPolyPMets.tsv');
genesInfo = fullfile(dataDir,'VolPolyPGenes.tsv');
rxnsCoeffs = fullfile(dataDir,'VolPolyPRxnsCoeffs.tsv');
rxnsInfo = fullfile(dataDir,'VolPolyPRxns.tsv');

model = curateMetsRxnsGenes(model, metsInfo, genesInfo, rxnsCoeffs, rxnsInfo);


%% DO NOT CHANGE OR REMOVE THE CODE BELOW THIS LINE.
% Show some metrics:
cd ../modelTests
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

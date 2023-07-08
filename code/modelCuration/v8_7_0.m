% This scripts applies curations to be applied on yeast-GEM release 8.6.3, to
% get to yeast-GEM release 8.7.0.
% Indicate which Issue/PR are addressed. If multiple curations are performed
% before a new release is made, just add the required code to this script. If
% more extensive coding is required, you can write a separate (generic) function
% that can be kept in the /code/modelCuration folder. Otherwise, try to use
% existing functions whenever possible. In particular /code/curateMetsRxnsGenes
% can do many types of curation.

%% Load yeast-GEM 8.6.3 (requires local yeast-GEM git repository)
cd ..
model = getEarlierModelVersion('8.6.3');
model.id='yeastGEM_develop';
dataDir=fullfile(pwd(),'..','data','modelCuration','v8_7_0');
cd modelCuration

%% Add new reactions based on KEGG and MetaCyc-derived reconstruction (PR #304)
metsInfo = fullfile(dataDir,'DBnewRxnsMets.tsv');
rxnsCoeffs = fullfile(dataDir,'DBnewRxnsCoeffs.tsv');
rxnsInfo = fullfile(dataDir,'DBnewRxnsRxns.tsv');
genesInfo = fullfile(dataDir,'DBnewRxnsGenes.tsv');
model = curateMetsRxnsGenes(model,metsInfo,genesInfo,rxnsCoeffs,rxnsInfo);

checkModelStruct(model,true,false)

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

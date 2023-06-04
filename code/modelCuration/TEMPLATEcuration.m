%% TEMPLATE FOR SCRIPT THAT CONSOLIDATES CURATIONS THAT ARE MADE TO A MODEL
% RELEASE. EACH RELEASE (e.g. yeast-GEM v8.6.1) WILL HAVE ITS OWN SCRIPT WHERE
% CURATIONS TO THAT MODEL RELEASE ARE CONSOLIDATED, MATCHING A DEDICATED FOLDER
% WITH DATA IN code/modelCuration/. SEE EARLIER SCRIPTS AS EXAMPLE. AFTER A NEW
% RELEASE IS MADE ON GITHUB, COPY-PASTE THIS TEMPLATE, REMOVE THIS LEADING
% CAPITALIZED TEXT, AND REPLACE $VERSION WITH THE VERSION NUMBER OF THE RECENT
% YEAST-GEM RELEASE (e.g. $VERSION --> 8.6.1).

% This scripts applies curations to be applied on yeast-GEM release $VERSION.
% Indicate which Issue/PR are addressed. If multiple curations are performed
% before a new release is made, just add the required code to this script. If
% more extensive coding is required, you can write a separate (generic) function
% that can be kept in the /code/modelCuration folder. Otherwise, try to use
% existing functions whenever possible. In particular /code/curateMetsRxnsGenes
% can do many types of curation.

%% Load yeast-GEM $VERSION (requires local yeast-GEM git repository)
cd ..
model = getEarlierModelVersion('$VERSION');
model.id='yeastGEM_develop';
dataDir=fullfile(pwd(),'..','data','modelCuration','v8.6.1'); %Modify to the correct version
cd modelCuration

%% Brief description of curation to be performed (PR #xxx) [include correct PR or Issue number]
% Potential longer description of curation to be performed.
% If any data files need to be loaded, keep these in the dedicated folder at
% data/modelCuration/v$VERSION

% Example [DELETE WHEN NOT USED]:

%% Curate gene association for transport rxns (PR #306)
% Add new genes
model       = curateMetsRxnsGenes(model,'none','../../data/modelCuration/curationsOnV8_6_0/transRxnNewGPRGenes.tsv');

% Change GPR relations
fid           = fopen('../../data/modelCuration/curationsOnV8_6_0/TransRxnNewGPR.tsv');
changegpr     = textscan(fid,'%q %q %q %q %q %q %q %q','Delimiter','\t','HeaderLines',1);
newGPR.ID     = changegpr{1};
newGPR.GPR    = changegpr{2};
fclose(fid);

model=changeGrRules(model,newGPR.ID,newGPR.GPR);

% Delete unused genes (if any)
model = deleteUnusedGenes(model);

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

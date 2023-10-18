% This scripts applies curations to be applied on yeast-GEM release 8.7.1, to
% get to yeast-GEM release 8.7.2.
% Indicate which Issue/PR are addressed. If multiple curations are performed
% before a new release is made, just add the required code to this script. If
% more extensive coding is required, you can write a separate (generic) function
% that can be kept in the /code/modelCuration folder. Otherwise, try to use
% existing functions whenever possible. In particular /code/curateMetsRxnsGenes
% can do many types of curation.

%% Load yeast-GEM 8.7.1 (requires local yeast-GEM git repository)
cd ..
codeDir=pwd();
model = getEarlierModelVersion('8.7.1');
model.id='yeastGEM_develop';
dataDir=fullfile(pwd(),'..','data','modelCuration','v8_7_2');
cd modelCuration

%% Correct dolichol-containing metFormulas
% While dolichol can have any number of isoprenoid units, in yeast-GEM it is
% defined as 4 units. This means that there is no need to keep R-subgroups as
% part of dolichol-derived metabolites to indicate that unspecified length.
% Define which metabolites are dolichol-derived and remove the R from their
% metabolite formula.
dolMets = getIndexes(model,{'s_3765','s_3767','s_3888','s_3911'},'mets');
model.metFormulas(dolMets) = regexprep(model.metFormulas(dolMets),'R','');

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

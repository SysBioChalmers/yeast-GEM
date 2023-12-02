% This scripts applies curations to be applied on yeast-GEM release 8.7.1, to
% get to yeast-GEM release 9.0.0.
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
model.version='';
%dataDir=fullfile(pwd(),'..','data','modelCuration','v9.0.0'); % No dataDir required for these curations
cd modelCuration

%% Correct inbalanced reactions, based on metFormulas
% While dolichol can have any number of isoprenoid units, in yeast-GEM it is
% defined as 4 units. This means that there is no need to keep R-subgroups as
% part of dolichol-derived metabolites to indicate that unspecified length.
% Define which metabolites are dolichol-derived and remove the R from their
% metabolite formula.
dolMets = getIndexes(model,{'s_3765','s_3767','s_3888','s_3911'},'mets');
model.metFormulas(dolMets) = regexprep(model.metFormulas(dolMets),'R','');

% r_4722 (polyphosphate polymerase) is unbalanced, 2 ADP is missing as product
model = changeRxns(model,'r_4722','2 ATP[c] + H2O[c] => H+[c] + polyphosphate[v] + 2 ADP[c]',3);

% r_4240 is unbalanced, but also has a generic reactant (protein asparagine)
% and describes a process (protein modification) that is out-of-scope of a
% metabolic model. Remove the reactions to resolve all 3 issues at once.
model = removeReactions(model,'r_4240',true,true,true);

% Some glycan metFormulas were incorrect absent. They were manually curated
% by drawing out the structures.
glycMets = getIndexes(model,{'s_3932','s_4003','s_4002'},'mets');
model.metFormulas(glycMets) = {'C50H82N4O37R','C38H70N2O36P2R2','C32H60N2O31P2R2'};

% r_0774 and r_0775 are unbalanced due to a missing H2O.
model = changeRxns(model,{'r_0774','r_0775'},...
    {'ATP[c] + H+[c] + nicotinate[c] + PRPP[c] => ADP[c] + diphosphate[c] + nicotinic acid D-ribonucleotide[c] + phosphate[c]',...
     'ATP[m] + H+[m] + nicotinate[m] + PRPP[m] => ADP[m] + diphosphate[m] + nicotinic acid D-ribonucleotide[m] + phosphate[m]'},3);

% r_4196 (NADH:ferricytochrome-b5 oxidoreductase) is unbalanced, NAD is missing as product
model = changeRxns(model,'r_4196','NADH[erm] + 2 Ferricytochrome b5[erm] <=> H+[erm] + 2 Ferrocytochrome b5[erm] + NAD[erm]',3);

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

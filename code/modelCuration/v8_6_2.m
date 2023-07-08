% This scripts applies curations to be applied on yeast-GEM release 8.6.1, to
% get to yeast-GEM release 8.6.2.
% Indicate which Issue/PR are addressed. If multiple curations are performed
% before a new release is made, just add the required code to this script. If
% more extensive coding is required, you can write a separate (generic) function
% that can be kept in the /code/modelCuration folder. Otherwise, try to use
% existing functions whenever possible. In particular /code/curateMetsRxnsGenes
% can do many types of curation.

%% Load yeast-GEM 8.6.1 (requires local yeast-GEM git repository)
cd ..
model = getEarlierModelVersion('8.6.1');
model.id='yeastGEM_develop';
dataDir=fullfile(pwd(),'..','data','modelCuration','v8_6_2');
cd modelCuration

%% Correct ATP synthase mitochondrial complex gene associations (PR #323)
model = changeGrRules(model, 'r_0226', ['Q0080 and Q0085 and Q0130 and ' ...
        'YBL099W and YBR039W and YDL004W and YDR298C and YDR377W and YJR121W ' ...
        'and YKL016C and YLR295C and YML081C-A and YPL078C and YPL271W and ' ...
        'YDR322C-A and YPR020W and YOL077W-A'],true);
model = deleteUnusedGenes(model);
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

% This scripts applies curations to be applied to generate yeast-GEM release 8.6.1.
% Indicate which Issue/PR are addressed

%% Load yeast-GEM 8.6.0 (requires local yeast-GEM git repository)
cd ..
model = getEarlierModelVersion('8.6.0');
model.id='yeastGEM_develop';

%% Curate complex annotation (PR #305)
% Add new genes
cd modelCuration
model       = curateMetsRxnsGenes(model,'none','../../data/modelCuration/v8_6_1/complexAnnotationGenes.tsv');

% Add gene standard name for new genes
fid = fopen('../../data/modelCuration/v8_6_1/complexAnnotation.tsv');
complexAnnot = textscan(fid,'%q %q %q %q %q %q %q','Delimiter','\t','HeaderLines',1);
fclose(fid);
newGPR.ID     = complexAnnot{1};
newGPR.GPR    = complexAnnot{3};
model=changeGrRules(model,newGPR.ID,newGPR.GPR);

% Delete unused genes (if any)
model = deleteUnusedGenes(model);

%% Curate gene association for transport rxns (PR #306)
% Add new genes
model       = curateMetsRxnsGenes(model,'none','../../data/modelCuration/v8_6_1/transRxnNewGPRGenes.tsv');

% Change GPR relations
fid           = fopen('../../data/modelCuration/v8_6_1/TransRxnNewGPR.tsv');
changegpr     = textscan(fid,'%q %q %q %q %q %q %q %q','Delimiter','\t','HeaderLines',1);
newGPR.ID     = changegpr{1};
newGPR.GPR    = changegpr{2};
fclose(fid);

model=changeGrRules(model,newGPR.ID,newGPR.GPR);

% Delete unused genes (if any)
model = deleteUnusedGenes(model);

%% Add new gene associations from databases (PR #313)
% Add new genes
model       = curateMetsRxnsGenes(model,'none','../../data/modelCuration/v8_6_1/newGPRsfromDBsGenes.tsv');

% Change GPR relations
fid           = fopen('../../data/modelCuration/v8_6_1/newGPRsfromDBs.tsv');
changegpr     = textscan(fid,'%q %q %q %q %q %q %q %q','Delimiter','\t','HeaderLines',1);
newGPR.ID     = changegpr{1};
newGPR.GPR    = changegpr{3};
fclose(fid);

model = changeGrRules(model,newGPR.ID,newGPR.GPR);

% Finding putative gene associations highlighted duplicated reactions:
% r_4566 is duplicate of r_4232, just in reverse direction
model = setParam(model,'lb','r_4232',-1000);
model = setParam(model,'rev','r_4232',1);
model = removeReactions(model,'r_4566',true,true,true);

% The following are duplicates, just in different compartment and with some
% different names for same metabolites (reactions were added based on
% metabolomics data, while isoleucine & valine biosynthesis is localized in
% the mitochondrion)
% r_4576 -> half-reaction of r_0096
% r_4577 -> duplicate of r_0352
% r_4578 -> duplicate of r_0096
% r_4579 -> duplicate of r_0097
% r_4580 -> duplicate of r_0096
model = removeReactions(model,...
    {'r_4576','r_4577','r_4578','r_4579','r_4580'},true,true,true);

% Delete unused genes (if any)
model = deleteUnusedGenes(model);

%% Define unique subsystems (Issue #11, PR #307)
fid           = fopen('../../data/modelCuration/v8_6_1/uniqueSubsystems.tsv');
fileInput     = textscan(fid,'%q %q %q %q %q %q %q','Delimiter','\t','HeaderLines',1);
fclose(fid);
subsystem.rxn = fileInput{1};
subsystem.sub = fileInput{5};

[a,b] = ismember(subsystem.rxn,model.rxns);
%Remove non-matching reactions
b(~a)=[]; subsystem.sub(~a)='';
for i=1:numel(b)
    model.subSystems{b(i),1}=subsystem.sub(i);
end

%% DO NOT CHANGE OR REMOVE THE CODE BELOW THIS LINE.
% Show some metrics:
cd ../modelTests
% disp('Run gene essentiality analysis')
% [new.accuracy,new.tp,new.tn,new.fn,new.fp] = essentialGenes(model);
% fprintf('Genes in model: %d\n',numel(model.genes));
% fprintf('Gene essentiality accuracy: %.4f\n', new.accuracy);
% fprintf('Gene essentiality TP: %d\n', numel(new.tp));
% fprintf('Gene essentiality TN: %d\n', numel(new.tn));
% fprintf('Gene essentiality FP: %d\n', numel(new.fp));
% fprintf('Gene essentiality FN: %d\n', numel(new.fn));
% fprintf('\nRun growth analysis\n')
% R2=growth(model);
% fprintf('R2 of growth prediction: %.4f\n', R2);

% Save model:
cd ..
saveYeastModel(model)
cd modelCuration

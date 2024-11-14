% This scripts applies curations to be applied on yeast-GEM release 9.0.1, to
% get to yeast-GEM release 9.0.2.

%% Load yeast-GEM 9.0.1 (requires local yeast-GEM git repository)
cd ..
codeDir=pwd();
model = getEarlierModelVersion('9.0.1');
model.id='yeastGEM_develop';
model.version='';
%dataDir=fullfile(pwd(),'..','data','modelCuration','v9.0.0'); % No dataDir required for these curations
cd modelCuration

%% Issue #371. Implement auxotrophy-based curation, as proposed by Han
% et al., 2024 doi:10.1016/j.synbio.2024.07.006. Not all changes are
% implemented as proposed by the authors.

% To resolve AAT2 causing asparate auxotrophy
model                       = setParam(model,'eq','r_0217',0);
rxnIdx                      = getIndexes(model,'r_0217','rxns');
model.rxnNotes{rxnIdx}      = 'No significant mitochondrial flux';
model.rxnReferences{rxnIdx} = '10.1016/j.femsyr.2004.09.008';

% To resolve ALD2 and ALD3 causing pantothenic acid auxotrophy
model                       = changeGrRules(model,'r_0172','YMR170C or YMR169C',true);
rxnIdx                      = getIndexes(model,'r_0172','rxns');
model.rxnReferences{rxnIdx} = '10.1093/genetics/163.1.69';

% To resolve CHO2 causing choline auxotrophy
CHO2rxns                    = {'r_2488','r_2489','r_2490','r_2491','r_2492','r_2493','r_2494','r_2495'};
model                       = changeGrRules(model,CHO2rxns,'YGR157W',true);
rxnIdx                      = getIndexes(model,CHO2rxns,'rxns');
model.rxnReferences(rxnIdx) = {'10.1016/0005-2760(90)90145-N'};
% While the curation is correct, CHO2 and OPI3 are not strictly essential
% genes, choline auxotrophy is only in relation to the production of
% phosphocholines. If no cholines are available, the phospholipid
% composition would be different. But because of the fixed biomass
% composition in yeast-GEM, CHO2 and OPI3 act as essential genes. This
% causes false positives in the essentialGenes function, as the genes not
% truly essential.

% To resolve CYS3 causing cysteine auxotrophy
model                       = setParam(model,'eq','r_0312',0);
rxnIdx                      = getIndexes(model,'r_0312','rxns');
model.rxnNotes{rxnIdx}      = 'Regulatory role only';
model.rxnReferences{rxnIdx} = '10.1111/j.1574-6968.2003.tb11531.x';

% Han et al. suggest to block TUM1 reaction r_4703, to resolve CYS4 as
% causing cysteine auxotrophy. However, there is good evidence that this
% reaction is taking place (doi:10.1093/femsyr/fow100). Instead, make the
% reaction irreversible.
model                       = setParam(model,'lb','r_4703',0);

% To resolve ERG10 causing ergosterol auxotrophy
model = removeReactions(model,'r_0559',true,true); % ERG13, mitochondrial
% only, see e.g. https://www.uniprot.org/uniprotkb/P54839/entry

% To resolve GFA1 causing D-glucosamine auxotrophy: add chitin to
% carbohydrate pseudoreaction and adjust other reactions proportionally
rxn  = getIndexes(model,'r_4048','rxns');
met1 = getIndexes(model,'s_0509','mets');%chitin
met2 = getIndexes(model,'s_0001','mets');%(1->3)-beta-D-glucan
met3 = getIndexes(model,'s_0004','mets');%(1->6)-beta-D-glucan
met4 = getIndexes(model,'s_1107','mets');%mannan
met5 = getIndexes(model,'s_0773','mets');%glycogen
met6 = getIndexes(model,'s_1520','mets');%trehalose
model.S(met1, rxn) = -0.02361;
model.S(met2, rxn) = -0.73914;
model.S(met3, rxn) = -0.24696;
model.S(met4, rxn) = -0.70204;
model.S(met5, rxn) = -0.35689;
model.S(met6, rxn) = -0.13655;
% And change gene association
model = changeGrRules(model,'r_0477','YKL104C',true);

% To resolve GSH1 causing glutatione auxotrophy: add glutathione to the
% cofactor pseudoreaction
rxn = getIndexes(model,'r_4598','rxns');
met = getIndexes(model,'s_0750','mets');
model.S(met, rxn) = -1e-06;
% The above change will make GSH1 and GSH2 essential genes. This will be
% reported in essentialGenes as false positives, as their essentiality is
% debated (see for instance doi:10.1016/S1567-1356(02)00081-8). We here
% argue that absence of glutathione very strongly affects cell growth and
% stress tolerance.

% To resolve ERG13/20 causing ergosterol auxotrophy; and HEM1/2/3/4/12 
% causing heme auxotrophy. Heme transport and exchange should be present to
% allow uptake if supplemented
metsToAdd.metNames      = 'heme a';
metsToAdd.compartments  = 'e';
model                   = addMets(model,metsToAdd,true,'s_');
model                   = addTransport(model,'c','e','heme a',true,false,'r_');
model                   = addExchangeRxns(model,'out',model.mets(end));
model.rxns(end)         = generateNewIds(model,'rxns','r_',1);
model.rxnNames(end)     = {'heme a exchange'};
model.subSystems{end}   = {'Exchange reaction'};

% To resolve MET3 causing methionine auxotrophy
rxnIdx = getIndexes(model,{'r_0223','r_0224','r_1026'},'rxns');
model.eccodes(rxnIdx) = {'2.7.7.53','2.7.7.53','2.7.7.5'};
% Han et al. suggest to block APA1 reaction r_1026, but their cited
% reference (doi:10.1007/BF00264692) contains no reasoning for this.
% Instead, there is good evidence for the function of this enzyme:
% (doi:10.1128/jb.171.12.6437-6445.1989), although the reaction is assayed
% in the reverse direction. So instead, invert the reaction (and keep it
% irreversible), this has the same effect for MET3 as blocking it.
model.S(:,rxnIdx(3)) = - model.S(:,rxnIdx(3));

% To resolve MET13 causing methionine auxotrophy
model = changeGrRules(model,'r_0080','(YGL125W and YPL023C) or YGL125W',true);

% To resolve MET17 causing methionine auxotrophy
model = setParam(model,'lb','r_0815',0);

% To resolve THI4 causing thiamine auxotrophy and HOM2 causing methionine
% and threonine auxotrophy
model = removeReactions(model,'r_2070',true,true);
model = removeReactions(model,'r_2071',true,true);
rxnsToAdd.rxns      = generateNewIds(model,'rxns','r_',2);
rxnsToAdd.equations = {'s_3910 + s_0803 => s_0423 + s_0293 + s_0456 + s_0794', ...
                       's_1003 + s_1198 + s_0841 => s_1216 + s_3910 + 3 s_0803 + s_0794'};
rxnsToAdd.rxnNames  = {'HET-P synthase (thiazole synthase)','adenylated thiazole synthase'};
rxnsToAdd.lb        = [0,0];
rxnsToAdd.ub        = [1000,1000];
rxnsToAdd.eccodes   = {'','2.4.2.60'};
rxnsToAdd.subSystems= {'Thiamine metabolism','Thiamine metabolism'};
rxnsToAdd.grRules   = {'YGR144W','YGR144W'};
rxnsToAdd.rxnMiriams{1} = struct('name',{{'kegg.reaction';'kegg.pathway';'metanetx.reaction'}},...
                                'value',{{'R10685';'sce00730';'MNXR139808'}});
rxnsToAdd.rxnMiriams{2} = struct('name',{{'kegg.reaction';'kegg.pathway';'metanetx.reaction'}},...
                                'value',{{'R10711';'sce00730';'MNXR139812'}});
rxnsToAdd.rxnReferences = {'10.1021/ja067606t','10.1021/ja067606t'};
rxnsToAdd.rxnConfidenceScores = [3,3];
model = addRxns(model,rxnsToAdd);

% To resolve URA2 causing uracil auxotrophy
model = changeGrRules(model,'r_0250','YJL130C or (YJR109C and YOR303W)',true);

%% Issue #372. Correct UniProt identifier of YCR024C and YCR024C-A
geneIdx = getIndexes(model,'YCR024C','genes');
model.geneMiriams{geneIdx}.value = {'P25345'};

%% DO NOT CHANGE OR REMOVE THE CODE BELOW THIS LINE.
% Show some metrics:
cd(fullfile(codeDir,'modelTests'))
disp('Run gene essentiality analysis')
[new.accuracy,new.tp,new.tn,new.fn,new.fp] = essentialGenes(model,true);
fprintf('Genes in model: %d\n',numel(model.genes));
fprintf('Gene essentiality accuracy: %.4f\n', new.accuracy);
fprintf('True non-essential genes: %d\n', numel(new.tp));
fprintf('True essential genes: %d\n', numel(new.tn));
fprintf('False non-essential genes: %d\n', numel(new.fp));
fprintf('False essential genes: %d\n', numel(new.fn));
fprintf('\nRun growth analysis\n')
R2=growth(model,true);
fprintf('R2 of growth prediction: %.4f\n', R2);

% Save model:
cd ..
saveYeastModel(model)
cd modelCuration

clear;
close all;
load('yeast-GEM.mat')

% this script makes corrections that should be applied to *all* models, not
% just the anaerobic. Adjustments to the anaerobic model are found in the
% anaerobicModel function


%=========================================================================

% Create a pathway for anaplerosis. Check alternatives? this one might work
% on other diacids also? If considered important for "flux reasons", dig 
% into this transporter literature
% 'r_1264' 'succinate transport'   'phosphate[mitochondrion] + succinate[cytoplasm]  -> phosphate[cytoplasm] + succinate[mitochondrion] '
%model.lb(findRxnIDs(model,'r_1264'))=-1000; %DIC1. I am not sure its actually a phosphate/succinate antiport. It not clear from the article. It could also be specific for both phosphate and succinate

%=========================================================================


% look for all proton symport/antiport reactions and make sure that they 
% only enter the cell. 

%'s_0794' H+ cytosol
%'s_0796' H+ extracellular
symporterIDs = intersect(find(model.S(find(strcmp(model.mets,'s_0796')),:)),find(model.S(find(strcmp(model.mets,'s_0794')),:)),'stable');
for i = 1:length(symporterIDs)
    if model.S(find(strcmp(model.mets,'s_0796')),symporterIDs(i))<0
        model.lb(symporterIDs(i))=0;
    else
        model.ub(symporterIDs(i))=0;
    end

    if find(strcmp(model.rxns,'r_1258'))==symporterIDs(i) % ignore the sodium transporter, without it, the model doesn't work
        model.lb(symporterIDs(i))=-1000;
        model.ub(symporterIDs(i))=1000;
    end
end

%=========================================================================

% This section balances reactions and ensures that a correct molecular
% weight can be calculated for the biomass

% Set the charge of all biomass components to 0, should be applied to *all* models, not just anaerobic.
model.metCharges(strcmp(model.mets,'s_3717'))=0; % Protein
model.metCharges(strcmp(model.mets,'s_3718'))=0; % Carbohydrate
model.metCharges(strcmp(model.mets,'s_3719'))=0; % RNA
model.metCharges(strcmp(model.mets,'s_3720'))=0; % DNA
model.metCharges(strcmp(model.mets,'s_3746'))=0; % Lipid backbone
model.metCharges(strcmp(model.mets,'s_3747'))=0; % Lipid chain
model.metCharges(strcmp(model.mets,'s_4205'))=0; % Cofactor
model.metCharges(strcmp(model.mets,'s_4206'))=0; % Ion

% Make the charge of K and Na 1+, should be applied to *all* models, not just anaerobic.
model.metCharges(strcmp(model.mets,'s_1373'))=1;
model.metCharges(strcmp(model.mets,'s_1374'))=1;
model.metCharges(strcmp(model.mets,'s_3776'))=1;
model.metCharges(strcmp(model.mets,'s_1437'))=1;
model.metCharges(strcmp(model.mets,'s_1438'))=1;
model.metCharges(strcmp(model.mets,'s_3775'))=1;

% Balance the charge of all biomass component pseudo reactions by adding the required amount of H+
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_4047')) = -sum(model.S(:,strcmp(model.rxns,'r_4047')).*model.metCharges,'omitnan'); % Protein
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_4049')) = -sum(model.S(:,strcmp(model.rxns,'r_4049')).*model.metCharges,'omitnan'); % RNA
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_4050')) = -sum(model.S(:,strcmp(model.rxns,'r_4050')).*model.metCharges,'omitnan'); % DNA
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_4598')) = -sum(model.S(:,strcmp(model.rxns,'r_4598')).*model.metCharges,'omitnan'); % Cofactor
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_4599')) = -sum(model.S(:,strcmp(model.rxns,'r_4599')).*model.metCharges,'omitnan'); % Ion


% Special case for SLIME rxns
model.metCharges(find(contains(model.metNames,'chain')+contains(model.metNames,'backbone'))) = 0;

% Now, based on the charge balance, find all the reactions that are
% imbalanced, add or remove hydrogen as necessary


% Balance the charge of all imbalanced SLIME reactions by adding the required amount of H+, 
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_3975')) = -sum(model.S(:,strcmp(model.rxns,'r_3975')).*model.metCharges,'omitnan'); % Protein
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_3976')) = -sum(model.S(:,strcmp(model.rxns,'r_3976')).*model.metCharges,'omitnan'); % Protein
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_3977')) = -sum(model.S(:,strcmp(model.rxns,'r_3977')).*model.metCharges,'omitnan'); % Protein
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_3978')) = -sum(model.S(:,strcmp(model.rxns,'r_3978')).*model.metCharges,'omitnan'); % Protein
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_4076')) = -sum(model.S(:,strcmp(model.rxns,'r_4076')).*model.metCharges,'omitnan'); % Protein
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_4077')) = -sum(model.S(:,strcmp(model.rxns,'r_4077')).*model.metCharges,'omitnan'); % Protein
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_4078')) = -sum(model.S(:,strcmp(model.rxns,'r_4078')).*model.metCharges,'omitnan'); % Protein
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_4079')) = -sum(model.S(:,strcmp(model.rxns,'r_4079')).*model.metCharges,'omitnan'); % Protein


% Make the charge of HS (hydrogen sulfide) -1
model.metCharges(strcmp(model.mets,'s_0841'))=-1;
model.metCharges(strcmp(model.mets,'s_3906'))=-1;
model.metCharges(strcmp(model.mets,'s_4263'))=-1;

% Balance the reaction r_4702, 'L-cysteine:2-oxoglutarate aminotransferase'
% by adding a proton as reactant
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_4702'))=-1;

% Balance the reaction r_4703, 'L-cysteine:2-oxoglutarate aminotransferase'
% by adding a proton as product
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_4703'))=1;

% Balance the reactions 'r_0774' and 'r_0775', 'NAPRtase' by removing H+ consumption
% and adding a H2O as a reactant
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_0774'))=0;
model.S(find(strcmp(model.mets,'s_0803')),strcmp(model.rxns,'r_0774'))=-1;
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_0775'))=0;
model.S(find(strcmp(model.mets,'s_0803')),strcmp(model.rxns,'r_0775'))=-1;

%=========================================================================
% This section focuses on individual reactions that have the wrong
% reversibility/direction/cofactor or should be completley removed

% rename r_0227, it is the plasma membrane ATPase, not a cytosolic ATPase
model.rxnNames(strcmp(model.rxns,'r_0227')) = {'ATPase, plasma membrane'};

% make sure both formate-THF ligases are reversible.
model.lb(strcmp(model.rxns,'r_0446')) = -1000; 

% ADE3 and MIS1, methylenetetrahydrofolate dehydrogenase (NADP+)  [EC 1.5.1.5]
% make irrevrersible
model.ub(strcmp(model.rxns,'r_0732')) = 0;
model.ub(strcmp(model.rxns,'r_0733')) = 0;

% There is no evidence for this PFK1 side reaction in yeast. Consider
% removing the reaction completley
model.ub(strcmp(model.rxns,'r_0887')) = 0; 

% TYR1 incorrectly annotated as using NAD, should be NADP 
model.S(find(strcmp(model.mets,'s_1212')),strcmp(model.rxns,'r_0939'))=0; %NADPH
model.S(find(strcmp(model.mets,'s_1207')),strcmp(model.rxns,'r_0939'))=0; %NADP
model.S(find(strcmp(model.mets,'s_1203')),strcmp(model.rxns,'r_0939'))=1; %NADH
model.S(find(strcmp(model.mets,'s_1198')),strcmp(model.rxns,'r_0939'))=-1;%NAD

% Make esterification reactions irreversible. positive deltaG
model.ub(strcmp(model.rxns,'r_4713')) = 0; %diethyl succinate
model.ub(strcmp(model.rxns,'r_4714')) = 0; %monoethyl succinate

% Make polyphosphate hydrolase and diphosphate transport over cell membrane
% both irreversible
model = setParam(model,'lb',{'r_4723','r_4724','r_4725'},0);
model = setParam(model,'lb','r_4460',0);

%=========================================================================


% Convert to anaerobic
model = anaerobicModel(model);

% Glucose uptake rate
model = setParam(model,'eq','r_1714',-23);

% Solve
res=solveLP(model,1);
FLUX=res.x;

v_glc=FLUX(getIndexes(model,'r_1714','rxns'),:);
v_eth=FLUX(getIndexes(model,'r_1761','rxns'),:);
v_CO2=FLUX(getIndexes(model,'r_1672','rxns'),:);
v_gly=FLUX(getIndexes(model,'r_1808','rxns'),:);
v_growth=FLUX(getIndexes(model,'r_4041','rxns'),:);
%%
%glycerol ethanol Co2
%4.5 ± 0.4  31 ± 2  38 ± 10
data=[4.5 31 38 0.36];
sim=[v_gly v_eth v_CO2 v_growth];
errorVal=[0.4 2 10 0.02];
b1=bar(data./data,'FaceAlpha',0.5);hold on;b2=bar(sim./data,'FaceAlpha',0.5);

hold on

er = errorbar([1 2 3 4],data./data,errorVal./data,errorVal./data);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  


legend({'data','simulation'});
ylabel('Relative value');
xticklabels({'Glycerol','Ethanol','CO2','Biomass'})

temp_model=model;
temp_model.metNames=strcat( model.metNames, repmat('[',length(model.mets),1),model.compNames(temp_model.metComps),repmat(']',length(model.mets),1) ); 


% Pack everything into a nice table
rxns_reacs=printRxnFormula(temp_model,'rxnAbbrList',model.rxns,'metNameFlag',1,'printFlag',0);
%[massImbalance, imBalancedMass, imBalancedCharge, imBalancedRxnBool, elements, missingFormulaeBool, balancedMetBool] = checkMassChargeBalance(model);
%tab=table(model.rxns,model.rxnNames,rxns_reacs,abs(FLUX./v_glc),FLUX,model.grRules,imBalancedRxnBool,imBalancedCharge,imBalancedMass);
tab=table(model.rxns,model.rxnNames,rxns_reacs,abs(FLUX./v_glc),FLUX,model.grRules);
tab = sortrows(tab,"Var4","descend");

printFluxes(model,FLUX,true)





%% Calculate formula and degree of reduciton of biomass

[mwRange,metFormulae,elements,metEle]=computeMetFormulae(model,'metMwRange','s_0450','fillMets','none','printLevel',0);
mwRange

% Print out the results of the 
Biomass_index = find(strcmp(model.metNames,'biomass'));

%Degree of reduction per element. Order of the elements 'C', 'H', 'O', 'N'
DR_per_ele = [4, 1, -2, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
DR_per_Cmol=sum(metEle(Biomass_index,:).*DR_per_ele)/metEle(Biomass_index,1)

Biomass_formula_Cmol=metEle(Biomass_index,:)/metEle(Biomass_index,1);

MW_per_Cmol_min = mwRange/metEle(Biomass_index,1)
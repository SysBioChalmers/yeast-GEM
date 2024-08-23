clear; close all
% Load the model and apply the corrections for *all* models
cd('../modelCuration/')
% run v_9.0.2.m;
v9_0_2;

% res=essentialGenes(model);
% res

%% Run growth tests
R2 = growth(model);

NLIM=0;
%% Convert to anaerobic
cd('../otherChanges/')
model = anaerobicModel(model,0);

cd('../modelTests/');
%% flux predictions
R2=anaerobic_flux_predictions(model);




%% Set glucose uptake rate and solve pFBA
model = setParam(model,'eq','r_1714',-23);
res=solveLP(model,1);
FLUX=res.x;


%% Retrieve data for the main products
v_glc=FLUX(getIndexes(model,'r_1714','rxns'),:);
v_eth=FLUX(getIndexes(model,'r_1761','rxns'),:);
v_CO2=FLUX(getIndexes(model,'r_1672','rxns'),:);
v_gly=FLUX(getIndexes(model,'r_1808','rxns'),:);
v_growth=FLUX(getIndexes(model,'r_4041','rxns'),:);


%% Show relative accuracy of main extracellular products
figure;
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

%% Pack flux results into table
temp_model=model;
temp_model.metNames=strcat( model.metNames, repmat('[',length(model.mets),1),model.compNames(temp_model.metComps),repmat(']',length(model.mets),1) ); 
rxns_reacs=printRxnFormula(temp_model,'rxnAbbrList',model.rxns,'metNameFlag',1,'printFlag',0);
tab=table(model.rxns,model.rxnNames,rxns_reacs,abs(FLUX./v_glc),FLUX,model.grRules);
tab = sortrows(tab,"Var4","descend");


%% Calculate formula and degree of reduciton of biomass
[mwRange,metFormulae,elements,metEle]=computeMetFormulae(model,'metMwRange','s_0450','fillMets','none','printLevel',0);
Biomass_index = find(strcmp(model.metNames,'biomass'));
%Degree of reduction per element. Order of the elements 'C', 'H', 'O', 'N'
DR_per_ele = [4, 1, -2, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
DR_per_Cmol=sum(metEle(Biomass_index,:).*DR_per_ele)/metEle(Biomass_index,1)
Biomass_formula_Cmol=metEle(Biomass_index,:)/metEle(Biomass_index,1);
MW_per_Cmol_min = mwRange/metEle(Biomass_index,1)




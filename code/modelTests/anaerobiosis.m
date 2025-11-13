clear; close all
% Load the model and apply the corrections for *all* models
%model = loadYeastModel;
cd('../modelCuration/')
v9_1_0;

%% Run growth tests
R2 = growth(model);

%% Convert to anaerobic
cd('../otherChanges/')
modelAn = anaerobicModel(model);

cd('../modelTests/');
%% flux predictions
R2=anaerobic_flux_predictions(modelAn);

%% Plot
plotAnaerobic(modelAn)
%% Set glucose uptake rate and solve pFBA
modelAn = setParam(modelAn,'eq','r_1714',-23);
res=solveLP(modelAn,1);
FLUX = res.x;
v_AStr = res.x(getIndexes(modelAn,'r_1115','rxns'));
v_ATPase = res.x(getIndexes(modelAn,'r_0227','rxns'));
v_glc = res.x(getIndexes(modelAn,'r_1714','rxns'));

%% Pack flux results into table
temp_model=modelAn;
rxns_reacs=constructEquations(temp_model,modelAn.rxns);
tab=table(modelAn.rxns,modelAn.rxnNames,rxns_reacs,abs(FLUX./v_glc),FLUX,modelAn.grRules);

[massImbalance, imBalancedMass, imBalancedCharge, imBalancedRxnBool, elements, missingFormulaeBool, balancedMetBool] = checkMassChargeBalance(modelAn);

tabImbalance = [tab,table(imBalancedRxnBool,imBalancedCharge,imBalancedMass)];
tabImbalance(~imBalancedRxnBool,:) = [];
% Filter out exchange and SLIME reactions, these are known to be unbalanced
tabImbalance(contains(tabImbalance.Var2,'exchange'),:) = [];
tabImbalance(contains(tabImbalance.Var2,'SLIME'),:) = [];


%% Calculate formula and degree of reduciton of biomass
% [mwRange,metFormulae,elements,metEle]=computeMetFormulae(modelAn,'metMwRange','s_0450','fillMets','none','printLevel',0);
% Biomass_index = find(strcmp(modelAn.metNames,'biomass'));
% %Degree of reduction per element. Order of the elements 'C', 'H', 'O', 'N'
% DR_per_ele = [4, 1, -2, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
% DR_per_Cmol=sum(metEle(Biomass_index,:).*DR_per_ele)/metEle(Biomass_index,1)
% Biomass_formula_Cmol=metEle(Biomass_index,:)/metEle(Biomass_index,1);
% MW_per_Cmol_min = mwRange/metEle(Biomass_index,1)
% 
% 


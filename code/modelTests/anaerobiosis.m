clear; close all
% Load the old and new model
model902 = getEarlierModelVersion('9.0.2');
model910 = loadYeastModel();

%% Run growth tests
funcs   = {@growth, @growthOld, @growth, @growthOld};
models  = {model910, model910, model902, model902};
titles  = {'v9.1.0 model, new anaerobic script', ...
           'v9.1.0 model, old anaerobic script', ...
           'v9.0.2 model, new anaerobic script', ...
           'v9.0.2 model, old anaerobic script'};

fig1 = figure('Name', 'Growth Comparison', 'Position', [100 100 1000 800]);
for k = 1:4
    subplot(2,2,k);
    funcs{k}(models{k});
    title(titles{k});
end
sgtitle('Growth: v9.1.0 vs v9.0.2, new vs old anaerobic script');
saveas(fig1, '..\..\data\testResults\v910_growth.png');

% R2 of growth prediction increased from 0.8256 to 0.9085 as result of all
% curations.

%% Convert to anaerobic
cd('../otherChanges/')
modelAn910 = anaerobicModel(model910);
modelAnOld910 = anaerobicModelOld(model910);
modelAn902 = anaerobicModel(model902);
modelAnOld902 = anaerobicModelOld(model902);

%% Anaerobic flux predictions
cd('../modelTests/');
models = {modelAn910, modelAnOld910, modelAn902, modelAnOld902};
titles  = {'v9.1.0 model, new anaerobic script', ...
           'v9.1.0 model, old anaerobic script', ...
           'v9.0.2 model, new anaerobic script', ...
           'v9.0.2 model, old anaerobic script'};
fig2 = figure('Name', 'Anaerobic Comparison', 'Position', [100 100 1000 800]);
for k = 1:4
    subplot(2,2,k);
    anaerobic_flux_predictions(models{k});
    title(titles{k});
end
sgtitle('Anaerobic flux predictions: v9.1.0 vs v9.0.2, new vs old anaerobic script');
saveas(fig2, '..\..\data\testResults\v910_anaerobic_fluxes.png');

% R2 of flux predictions increased from 0.7858 to 0.9175 as a response to
% all curations

%% Plot
models = {modelAn910, modelAnOld910, modelAn902, modelAnOld902};
titles  = {'v9.1.0 model, new anaerobic script', ...
           'v9.1.0 model, old anaerobic script', ...
           'v9.0.2 model, new anaerobic script', ...
           'v9.0.2 model, old anaerobic script'};
fig3 = figure('Name', 'Anaerobic Comparison', 'Position', [100 100 1000 800]);
for k = 1:4
    subplot(2,2,k);
    plotAnaerobic(models{k});
    title(titles{k});
end
sgtitle('Anaerobic fluxes (Sjöberg data): v9.1.0 vs v9.0.2, new vs old anaerobic script');
saveas(fig3, '..\..\data\testResults\v910_anaerobic_sjoberg.png');

%% Set glucose uptake rate and solve pFBA
temp_model = setParam(modelAn910,'eq','r_1714',-23);
res=solveLP(temp_model,1); FLUX = res.x;

v_AStr = res.x(getIndexes(temp_model,'r_1115','rxns')); % Ammonium exchange
v_ATPase = res.x(getIndexes(temp_model,'r_0227','rxns')); % ATPase
v_glc = res.x(getIndexes(temp_model,'r_1714','rxns')); % Glucose uptake
[v_AStr, v_ATPase, v_glc]

%% Pack flux results into table
rxns_reacs=constructEquations(temp_model,temp_model.rxns);
tab=table(temp_model.rxns,temp_model.rxnNames,rxns_reacs,abs(FLUX./v_glc),FLUX,temp_model.grRules);

[massImbalance, imBalancedMass, imBalancedCharge, imBalancedRxnBool, elements, missingFormulaeBool, balancedMetBool] = checkMassChargeBalance(temp_model);

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


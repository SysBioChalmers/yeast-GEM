%This function is to  check whether adding new reaction will lead to infinite
%ATP and NADH production
%the input is model and new reaction ID lists
%output is a result with pass or error


%
% Feiran Li
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [EnergyResults,RedoxResults] = CheckEnergyProduction(model,rxn,EnergyResults,RedoxResults)
[~,rxnID] = ismember(rxn,model.rxns);
%EnergyResults = {};
%RedoxResults = {};
model_test = model;
model_test = minimal_Y6(model_test);
%Add/change ATP production reaction:
%            ATP    +    H2O    ->  ADP     +   H+      +  PO4
mets  = {'s_0434[c]','s_0803[c]','s_0394[c]','s_0794[c]','s_1322[c]'};
coefs = [-1,-1,1,1,1];
model_test = addReaction(model_test,{'GenerateATP','leaktest1'}, ...
                    mets,coefs,false,0,1000);
model_test = changeObjective(model_test, model_test.rxns(end), 1);
sol = optimizeCbModel(model_test);
if sol.obj <= 360 && sol.obj > 0 %later can be changed to the experimental value
  EnergyResults = [EnergyResults; model.rxns(rxnID),'pass',num2str(sol.obj)];
elseif sol.obj > 360
  EnergyResults = [EnergyResults; model.rxns(rxnID),'Fail',num2str(sol.obj)];
else
  EnergyResults = [EnergyResults; model.rxns(rxnID),'error','error'];
end


model_test = model;
model_test = minimal_Y6(model_test);
%Add/change NADH production reaction:
                    %            NADH[c] + h[c] =>  NAD[c] 
                    mets  = {'s_1203[c]','s_0794[c]','s_1198[c]'};
                    coefs = [-1,-1,1];
                    model_test = addReaction(model_test,{'GenerateNADH','leaktest12'}, ...
                                        mets,coefs,false,0,1000);
model_test = changeObjective(model_test, model_test.rxns(end), 1);
sol = optimizeCbModel(model_test);
if sol.obj <= 120 && sol.obj > 0 %later can be changed to the experimental value
    RedoxResults = [RedoxResults; model.rxns(rxnID),'pass',num2str(sol.obj)];
elseif sol.obj > 120
    RedoxResults = [RedoxResults; model.rxns(rxnID),'Fail',num2str(sol.obj)];
elseif sol.obj <= 0
    RedoxResults = [RedoxResults; model.rxns(rxnID),'error','error'];
end





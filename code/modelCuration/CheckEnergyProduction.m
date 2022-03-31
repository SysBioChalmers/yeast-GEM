% This function is to  check whether adding new reaction will lead to infinite
% ATP and NADH production
% the input is model and new reaction ID lists
% output is a result with pass or error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [energyResults,redoxResults] = CheckEnergyProduction(model,rxn,energyResults,redoxResults)

[~,rxnID] = ismember(rxn,model.rxns);

if rxnID ~= 0
    model_test = model;
    model_test = minimal_Y6(model_test);
    %Add/change ATP production reaction:
    %            ATP    +    H2O    ->  ADP     +   H+      +  PO4
    mets  = {'s_0434[c]','s_0803[c]','s_0394[c]','s_0794[c]','s_1322[c]'};
    coefs = [-1,-1,1,1,1];
    model_test = addReaction(model_test,{'GenerateATP','leaktest1'}, ...
                             mets,coefs,false,0,1000);
    model_test = changeObjective(model_test,'GenerateATP', 1);
    sol = optimizeCbModel(model_test);
    if sol.f <= 360 && sol.f > 0 %later can be changed to the experimental value
        energyResults = [energyResults; model.rxns(rxnID),'pass',num2str(sol.f)];
    elseif sol.f > 360
        energyResults = [energyResults; model.rxns(rxnID),'Fail',num2str(sol.f)];
    else
        energyResults = [energyResults; model.rxns(rxnID),'error','error'];
    end
    
    model_test = model;
    model_test = minimal_Y6(model_test);
    %Add/change NADH production reaction:
    %            NADH[c] + H[c] =>  NAD[c]
    mets  = {'s_1203[c]','s_0794[c]','s_1198[c]'};
    coefs = [-1,-1,1];
    model_test = addReaction(model_test,{'GenerateNADH','leaktest12'}, ...
                             mets,coefs,false,0,1000);
    model_test = changeObjective(model_test, model_test.rxns(end), 1);
    sol = optimizeCbModel(model_test);
    if sol.f <= 120 && sol.f > 0 %later can be changed to the experimental value
        redoxResults = [redoxResults; model.rxns(rxnID),'pass',num2str(sol.obj)];
    elseif sol.obj > 120
        redoxResults = [redoxResults; model.rxns(rxnID),'Fail',num2str(sol.obj)];
    elseif sol.obj <= 0
        redoxResults = [redoxResults; model.rxns(rxnID),'error','error'];
    end
else
    redoxResults = [redoxResults; {'alreadlyexist'},'skip','skip'];
end

end

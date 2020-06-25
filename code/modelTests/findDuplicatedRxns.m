function findDuplicatedRxns(model)
% findDuplicatedRxns
%   Find and print reactions that have the same stoichiometry (forwards or
%   backwards).
%
%   Input:
%   model   genome-scale model
%
%   Usage: findDuplicatedRxns(model)
%

model_r = ravenCobraWrapper(model);
for i = 1:length(model.rxns)-1
    for j = i+1:length(model.rxns)
        if isequal(model.S(:,i),model.S(:,j)) || isequal(model.S(:,i),-model.S(:,j))
            printRxnFormula(model,model.rxns(i),true,true,true);
            disp(['Name: ' model.rxnNames{i} ' - GPR: ' model_r.grRules{i} ' - LB=' num2str(model.lb(i)) ' - UB=' num2str(model.ub(i))])
            printRxnFormula(model,model.rxns(j),true,true,true);
            disp(['Name: ' model.rxnNames{j} ' - GPR: ' model_r.grRules{j} ' - LB=' num2str(model.lb(j)) ' - UB=' num2str(model.ub(j))])
            disp(" ")
        end
    end
end

end

function findDuplicatedRxns(model)
% findDuplicatedRxns  yeast-GEM shim — delegates to RAVEN's
% findDuplicateRxns and prints each pair in the legacy format.
%
%   For every (i, j) pair of reactions sharing stoichiometry (in
%   either direction), prints two lines with name / GPR / lb / ub —
%   matching the pre-refactor output verbatim.
%
% Usage: findDuplicatedRxns(model)

pairs = findDuplicateRxns(model);
for k = 1:size(pairs, 1)
    i = pairs(k, 1);
    j = pairs(k, 2);
    constructEquations(model, model.rxns(i));
    disp(['Name: ' model.rxnNames{i} ' - GPR: ' model.grRules{i} ...
          ' - LB=' num2str(model.lb(i)) ' - UB=' num2str(model.ub(i))])
    constructEquations(model, model.rxns(j));
    disp(['Name: ' model.rxnNames{j} ' - GPR: ' model.grRules{j} ...
          ' - LB=' num2str(model.lb(j)) ' - UB=' num2str(model.ub(j))])
    disp(" ")
end
end

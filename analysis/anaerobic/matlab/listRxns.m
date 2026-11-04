function listRxns(names)
% listRxns  Print every v9.1.0 reaction that involves a metabolite whose name
%   matches one of the given names exactly.
paths = setupStudy('3.0.0-rc1');
m = loadTaggedModel(paths.yeast, '9.1.0');
for n = names
    idx = find(any(m.S(strcmp(m.metNames, n{1}), :) ~= 0, 1));
    eqs = constructEquations(m, m.rxns(idx));
    fprintf('\n== %s\n', n{1});
    for k = 1:numel(idx)
        i = idx(k);
        fprintf('%s [%g,%g] %s | %s | %s\n', m.rxns{i}, m.lb(i), m.ub(i), m.rxnNames{i}, m.grRules{i}, eqs{k});
    end
end
end

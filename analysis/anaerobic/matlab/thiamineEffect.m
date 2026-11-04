function thiamineEffect()
% thiamineEffect  v9.1.1 anaerobic (Sjoberg batch): growth, and anaerobic
%   single-gene deletions whose predicted viability changes when thiamine
%   uptake is opened and Thi5 (r_1603) is blocked.
paths = setupStudy('3.0.0-rc1');
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
a = variantModel(paths, 'v911_new');
b = setParam(setParam(a, 'ub', 'r_1603', 0), 'lb', 'r_2067', -1000);
gr = cell(1,2); k = 0;
for m = {a, b}
    k = k + 1;
    mm = setParam(setParam(m{1}, 'eq', 'r_1714', -23), 'obj', 'r_4041', 1);
    s = solveLP(mm, 1);
    fprintf('model %d: growth %.6f, r_1603 %.2e, thiamine uptake %.2e\n', k, s.x(strcmp(mm.rxns,'r_4041')), s.x(strcmp(mm.rxns,'r_1603')), -s.x(strcmp(mm.rxns,'r_2067')));
    [genes, ~, ~, ~, ratio] = findGeneDeletions(mm, 'sgd', 'fba');
    v = ones(numel(mm.genes),1); v(genes) = ratio; gr{k} = containers.Map(mm.genes, num2cell(v));
end
g = keys(gr{1});
for i = 1:numel(g)
    if abs(gr{1}(g{i}) - gr{2}(g{i})) > 1e-3
        fprintf('  %-10s without %.3f  with %.3f\n', g{i}, gr{1}(g{i}), gr{2}(g{i}));
    end
end
end

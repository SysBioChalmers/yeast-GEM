function gpdSinksGH()
% gpdSinksGH  gpd1 gpd2 in package G+H (GCV open, Hom6 NADPH only): which
%   reactions carry the redox that glycerol would carry, and how growth
%   responds when each candidate sink is closed in turn (bound-free checks:
%   exchange closed or reaction removed, one at a time and combined).

paths = setupStudy('3.0.0-rc1');
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
m = loadTaggedModel(paths.yeast, '9.1.0');
m = setParam(setParam(m, 'ub', {'r_0501','r_0507','r_0509'}, 1000), 'ub', 'r_0546', 0);
m = anaerobicModel(m);
m = removeGenes(m, {'YDL022W','YOL059W'});
m = setParam(setParam(m, 'eq', 'r_1714', -19.1), 'obj', 'r_4041', 1);
s = solveLP(m, 1);
fprintf('gpd1 gpd2, G+H: growth %.4f\n', s.x(strcmp(m.rxns,'r_4041')));
for n = {'succinate','fumarate','FADH2','propanol'}
    for i = find(strcmp(m.metNames, n{1}))'
        rx = find(m.S(i,:) ~= 0 & abs(s.x') > 1e-3);
        eqs = constructEquations(m, m.rxns(rx));
        for k = 1:numel(rx)
            fprintf('  [%s %s] %-8s %+8.3f %s | %s\n', n{1}, m.comps{m.metComps(i)}, m.rxns{rx(k)}, ...
                full(m.S(i,rx(k)))*s.x(rx(k)), m.rxnNames{rx(k)}, eqs{k});
        end
    end
end
closeSets = {'succinate', {'r_2056'}; 'propanol', {'r_4659'}; 'succinate+propanol', {'r_2056','r_4659'}};
for c = 1:size(closeSets,1)
    mm = setParam(m, 'ub', closeSets{c,2}, 0);
    sc = solveLP(mm, 1);
    if isempty(sc.x), fprintf('closed %s: infeasible\n', closeSets{c,1}); continue; end
    [~, exch] = getExchangeRxns(mm);
    big = exch(sc.x(exch) > 0.05 & ~ismember(mm.rxns(exch), {'r_1761','r_1672','r_2100','r_1832','r_4041'}));
    fprintf('closed %s: growth %.4f; exports %s\n', closeSets{c,1}, sc.x(strcmp(mm.rxns,'r_4041')), ...
        strjoin(strcat(mm.rxnNames(big)', '=', arrayfun(@(v) sprintf('%.2f',v), sc.x(big)', 'UniformOutput', false)), ', '));
end
end

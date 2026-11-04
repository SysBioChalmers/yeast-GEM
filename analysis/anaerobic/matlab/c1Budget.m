function c1Budget()
% c1Budget  Where the one-carbon units go in v9.1.0 + new constraints
%   (Sjoberg batch, pFBA): every reaction that consumes or produces
%   5,10-methylene-THF, 5,10-methenyl-THF, 10-formyl-THF, 5-methyl-THF or
%   S-adenosyl-L-methionine, with its net contribution, per compartment.

paths = setupStudy('3.0.0-rc1');
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
m = anaerobicModel(loadTaggedModel(paths.yeast, '9.1.0'));
m = setParam(m, 'eq', 'r_1714', -23);
m = setParam(m, 'obj', 'r_4041', 1);
sol = solveLP(m, 1);
c1 = {'5,10-methylenetetrahydrofolate','5,10-methenyl-THF','10-formyl-THF','5-methyltetrahydrofolate','S-adenosyl-L-methionine'};
fprintf('C1 carriers in the model: %s\n', strjoin(unique(m.metNames(contains(m.metNames, {'THF','tetrahydrofol','S-adenosyl'}))), ' | '));
for n = c1
    idx = find(strcmp(m.metNames, n{1}));
    if isempty(idx), fprintf('\n(no metabolite named %s)\n', n{1}); continue; end
    fprintf('\n=== %s\n', n{1});
    for i = idx'
        rx = find(m.S(i,:) ~= 0 & abs(sol.x') > 1e-6);
        [~, o] = sort(abs(full(m.S(i,rx)) .* sol.x(rx)'), 'descend');
        for j = rx(o)
            fprintf('  [%s] %-8s %+8.4f  %s\n', m.comps{m.metComps(i)}, m.rxns{j}, full(m.S(i,j))*sol.x(j), m.rxnNames{j}(1:min(60,end)));
        end
    end
end
end

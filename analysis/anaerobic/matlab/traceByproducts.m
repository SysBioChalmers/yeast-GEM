function traceByproducts()
% traceByproducts  Why does v9.1.0 + new constraints (Sjoberg batch, glucose
%   23) need to secrete glycine and glycolaldehyde?
%   1. Minimum secretion of each at optimal growth (forced or pFBA's choice),
%      and growth when the secretion is blocked.
%   2. In the pFBA solution, every reaction producing or consuming the
%      metabolite (all compartments), and the net production per reaction.
%   3. The same with the curation-blocked reactions (glycine cleavage
%      r_0501/r_0507/r_0509, glutamate synthase r_0472) reopened.

paths = setupStudy('3.0.0-rc1');
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
m = anaerobicModel(loadTaggedModel(paths.yeast, '9.1.0'));
m = setParam(m, 'eq', 'r_1714', -23);
m = setParam(m, 'obj', 'r_4041', 1);
sol = solveLP(m, 1);
mu = sol.x(strcmp(m.rxns, 'r_4041'));
fprintf('max growth %.5f\n', mu);

byproducts = {'r_1810', 'L-glycine'; 'r_1814', 'glycolaldehyde'};
for b = 1:size(byproducts,1)
    ex = byproducts{b,1};
    fprintf('\n=== %s (%s): pFBA secretion %.4f\n', byproducts{b,2}, ex, sol.x(strcmp(m.rxns, ex)));
    % Minimum secretion at (near) optimal growth.
    mm = setParam(m, 'lb', 'r_4041', 0.99999*mu);
    mm = setParam(mm, 'obj', ex, -1);                  % minimise secretion
    s = solveLP(mm);
    fprintf('  minimum secretion at optimal growth: %.4f\n', s.x(strcmp(mm.rxns, ex)));
    % Growth with the secretion blocked.
    mb = setParam(m, 'ub', ex, 0);
    s = solveLP(mb);
    if isempty(s.x), fprintf('  growth with secretion blocked: infeasible\n');
    else, fprintf('  growth with secretion blocked: %.5f (%.1f%% of max)\n', s.x(strcmp(mb.rxns,'r_4041')), 100*s.x(strcmp(mb.rxns,'r_4041'))/mu); end
    % Producers and consumers in the pFBA solution.
    mets = find(strcmp(m.metNames, byproducts{b,2}));
    fprintf('  reactions using %s (pFBA flux x coefficient, per compartment):\n', byproducts{b,2});
    for i = mets'
        rx = find(m.S(i,:) ~= 0 & abs(sol.x') > 1e-7);
        for j = rx
            eq = constructEquations(m, m.rxns(j));
            fprintf('    [%s] %-8s %+8.4f  %-40s %s\n', m.comps{m.metComps(i)}, m.rxns{j}, full(m.S(i,j))*sol.x(j), ...
                m.rxnNames{j}(1:min(40,end)), eq{1}(1:min(110,end)));
        end
    end
end

fprintf('\n=== With the curation-blocked reactions reopened\n');
blocked = {'r_0501','r_0507','r_0509','r_0472'};
for k = 1:numel(blocked)
    fprintf('  %s %s: [%g, %g] %s\n', blocked{k}, m.rxnNames{strcmp(m.rxns,blocked{k})}, ...
        m.lb(strcmp(m.rxns,blocked{k})), m.ub(strcmp(m.rxns,blocked{k})), strjoin(constructEquations(m, blocked(k)), ''));
end
for set = {{'r_0501','r_0507','r_0509'}, {'r_0472'}, {'r_0501','r_0507','r_0509','r_0472'}}
    mo = m;
    for r = set{1}, mo = setParam(mo, 'ub', r{1}, 1000); end
    s = solveLP(mo, 1);
    fprintf('  reopened %s: growth %.5f, glycine secretion %.4f, glycolaldehyde %.4f, glycerol %.4f\n', ...
        strjoin(set{1}, '+'), s.x(strcmp(mo.rxns,'r_4041')), s.x(strcmp(mo.rxns,'r_1810')), ...
        s.x(strcmp(mo.rxns,'r_1814')), s.x(strcmp(mo.rxns,'r_1808')));
end
end

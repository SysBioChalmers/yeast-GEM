function gcvLoop()
% gcvLoop  What does v9.1.0 + new anaerobic constraints do with glycine
%   cleavage (GCV) reopened, and which bound-free changes stop it?
%   1. Reactions whose pFBA flux changes most when GCV is reopened
%      (Sjoberg batch, glucose 23), with equations: the pathway that
%      replaces glycerol formation.
%   2. Bound-free variants, each scored on growth, glycerol, GCV flux and
%      glycine export for the wild type and gpd1 gpd2.

paths = setupStudy('3.0.0-rc1');
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
m0 = anaerobicModel(loadTaggedModel(paths.yeast, '9.1.0'));
m0 = setParam(setParam(m0, 'eq', 'r_1714', -23), 'obj', 'r_4041', 1);
gcv = {'r_0501','r_0507','r_0509'};
open = setParam(m0, 'ub', gcv, 1000);

s0 = solveLP(m0, 1); s1 = solveLP(open, 1);
d = s1.x - s0.x;
[~, o] = sort(abs(d), 'descend');
o = o(abs(d(o)) > 0.02);
eqs = constructEquations(m0, m0.rxns(o));
fprintf('Reactions changing when GCV is reopened (pFBA, glucose 23):\n');
fprintf('%-8s %8s %8s %8s | name | equation\n', 'rxn', 'blocked', 'open', 'diff');
for k = 1:numel(o)
    i = o(k);
    fprintf('%-8s %8.3f %8.3f %+8.3f | %s | %s\n', m0.rxns{i}, s0.x(i), s1.x(i), d(i), m0.rxnNames{i}, eqs{k});
end

v = {'released (GCV blocked)', m0};
v(end+1,:) = {'GCV open', open};
v(end+1,:) = {'GCV open, Shm1 serine->glycine only', setParam(open, 'lb', 'r_0503', 0)};
v(end+1,:) = {'GCV open, Shm1 glycine->serine only', setParam(open, 'ub', 'r_0503', 0)};
v(end+1,:) = {'GCV open, threonine aldolase blocked', setParam(open, 'ub', 'r_1040', 0)};
v(end+1,:) = {'GCV open, glycine export closed', setParam(open, 'ub', 'r_1810', 0)};
v(end+1,:) = {'GCV open, Shm1 ser->gly only, gly export closed', setParam(setParam(open, 'lb', 'r_0503', 0), 'ub', 'r_1810', 0)};
v(end+1,:) = {'GCV open, Hom6 NADPH only', setParam(open, 'ub', 'r_0546', 0)};
v(end+1,:) = {'released, Hom6 NADPH only', setParam(m0, 'ub', 'r_0546', 0)};

rx = {'r_4041','r_1808','r_1810','r_0501','r_0502','r_0503','r_1040','r_0546','r_4659'};
fprintf('\n%-50s %-9s %s\n', 'variant', 'strain', strjoin(rx, '   '));
for k = 1:size(v,1)
    for st = {'WT', {}; 'gpd1gpd2', {'YDL022W','YOL059W'}}'
        mm = v{k,2}; if ~isempty(st{2}), mm = removeGenes(mm, st{2}); end
        s = solveLP(mm, 1);
        if isempty(s.x), fprintf('%-50s %-9s infeasible\n', v{k,1}, st{1}); continue; end
        f = arrayfun(@(r) s.x(strcmp(mm.rxns, r{1})), rx);
        fprintf('%-50s %-9s %s\n', v{k,1}, st{1}, sprintf('%7.4f ', f));
    end
end
end

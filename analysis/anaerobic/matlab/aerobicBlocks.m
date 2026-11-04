function aerobicBlocks()
% aerobicBlocks  Growth benchmark (32 chemostats, v9.1.0, current growth.m set-up)
%   with each base-model choice of 9.1.0 undone one at a time and all together:
%   glutamate synthase r_0472 open, glycine cleavage r_0501/r_0507/r_0509 open,
%   Gcy1 r_0487 reversible, biomass reducing-equivalent term (3 NADPH) removed.
%   The choices sit in the base model file, so they act in aerobic and anaerobic
%   simulations alike. Prints growth R2 and the mean relative error per group.

paths = setupStudy('3.0.0-rc1');
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
base = loadTaggedModel(paths.yeast, '9.1.0');
gcv = {'r_0501','r_0507','r_0509'};
glt = @(m) setParam(m, 'ub', 'r_0472', 1000);
gcvo = @(m) setParam(m, 'ub', gcv, 1000);
gcy = @(m) setParam(m, 'lb', 'r_0487', -1000);
bio = @(m) removeBiomassTerm(m);
V = {'as released', base; 'Glt1 open', glt(base); 'GCV open', gcvo(base); 'Gcy1 reversible', gcy(base); ...
     'biomass NADPH term removed', bio(base); 'all four undone', bio(gcy(gcvo(glt(base))))};
grp = [repmat({'N-lim aerobic'},9,1); repmat({'C-lim aerobic'},11,1); repmat({'C-lim anaerobic'},6,1); repmat({'N-lim anaerobic'},6,1)];
names = unique(grp, 'stable');
fprintf('%-28s %8s %8s %12s %12s %12s %12s\n', 'variant', 'R2', 'MRE', names{:});
for k = 1:size(V,1)
    figure('Visible','off'); [~, e, g] = growthNewSim(V{k,2}); close all
    o = double(e(:,4)); p = g(:,4);
    r2 = 1 - sum((o-p).^2)/sum((o-mean(o)).^2); rel = abs(p-o)./o;
    fprintf('%-28s %8.4f %8.4f', V{k,1}, r2, mean(rel));
    for n = names', fprintf(' %12.4f', mean(rel(strcmp(grp, n{1})))); end
    fprintf('\n');
end
end

function m = removeBiomassTerm(m)
bio = strcmp(m.rxns, 'r_4041');
[~, mi] = ismember({'s_1212','s_1207','s_0794'}, m.mets);   % NADPH[c], NADP[c], H+[c]
m.S(mi, bio) = m.S(mi, bio) + [3; -3; 3];
assert(abs(m.S(mi(1), bio)) < 1e-9, 'NADPH not removed from r_4041');
end

function tyr1Variants()
% tyr1Variants  Prephenate dehydrogenase r_0939 (TYR1) with NAD/NADH (v9.1.0 as
%   released) and with NADP/NADPH (v9.0.2 form): Sjoberg batch (glucose 23),
%   anaerobic 13C flux benchmark, growth benchmark, aerobic growth on the
%   minimal medium of the yeast-GEM growth check.

paths = setupStudy('3.0.0-rc1');
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
base = loadTaggedModel(paths.yeast, '9.1.0');
r = strcmp(base.rxns, 'r_0939');
ix = @(id) find(strcmp(base.mets, id));
nadp = base;
nadp.S([ix('s_1198') ix('s_1203') ix('s_1207') ix('s_1212')], r) = [0; 0; -1; 1];   % NAD, NADH off; NADP -1, NADPH +1
assert(full(nadp.S(ix('s_1207'), r)) == -1 && full(nadp.S(ix('s_1198'), r)) == 0);
V = {'NAD/NADH (released)', base; 'NADP/NADPH', nadp};
grp = [repmat({'N-lim aerobic'},9,1); repmat({'C-lim aerobic'},11,1); repmat({'C-lim anaerobic'},6,1); repmat({'N-lim anaerobic'},6,1)];
names = unique(grp, 'stable');
fprintf('%-22s %7s %8s %8s | %7s %7s | %7s %7s', 'form', 'growth', 'glycerol', 'ethanol', 'fluxR2', 'fluxMRE', 'grR2', 'grMRE');
fprintf(' %12s', names{:}); fprintf('\n');
for k = 1:size(V,1)
    m = V{k,2}; an = anaerobicModel(m);
    an = setParam(setParam(an, 'eq', 'r_1714', -23), 'obj', 'r_4041', 1);
    s = solveLP(an, 1); x = @(q) s.x(strcmp(an.rxns, q));
    figure('Visible','off'); [r2f, mre] = fluxPredSim(anaerobicModel(m)); close all
    figure('Visible','off'); [~, e, g] = growthNewSim(m); close all
    o = double(e(:,4)); p = g(:,4); rel = abs(p-o)./o;
    fprintf('%-22s %7.4f %8.3f %8.2f | %7.4f %7.4f | %7.4f %7.4f', V{k,1}, x('r_4041'), x('r_1808'), x('r_1761'), r2f, mre, ...
        1 - sum((o-p).^2)/sum((o-mean(o)).^2), mean(rel));
    for n = names', fprintf(' %12.4f', mean(rel(strcmp(grp, n{1})))); end
    fprintf('\n');
end
end

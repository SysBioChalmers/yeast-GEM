function biomassTermVariants()
% biomassTermVariants  The 3 mmol/gDW reducing-equivalent term of the 9.1.0 biomass
%   reaction r_4041 in three forms: NADPH (as released), removed, NADH (as the
%   manuscript and Sjoberg 2024 describe it). Per form: Sjoberg batch (glucose 23)
%   growth, glycerol and ethanol; anaerobic 13C flux benchmark (R2, mean relative
%   error); growth benchmark (R2, mean relative error, per group).

paths = setupStudy('3.0.0-rc1');
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
base = loadTaggedModel(paths.yeast, '9.1.0');
bio = strcmp(base.rxns, 'r_4041');
ix = @(ids) cellfun(@(s) find(strcmp(base.mets, s)), ids);
nadph = ix({'s_1212','s_1207','s_0794'});     % NADPH[c], NADP[c], H+[c]
nadh  = ix({'s_1203','s_1198','s_0794'});     % NADH[c], NAD[c], H+[c]
assert(abs(base.S(nadph(1), bio) + 3) < 1e-9, 'released biomass has no -3 NADPH');
none = base; none.S(nadph, bio) = none.S(nadph, bio) + [3; -3; 3];
asNadh = none; asNadh.S(nadh, bio) = asNadh.S(nadh, bio) + [-3; 3; -3];
V = {'NADPH (released)', base; 'removed', none; 'NADH (as published)', asNadh};
grp = [repmat({'N-lim aerobic'},9,1); repmat({'C-lim aerobic'},11,1); repmat({'C-lim anaerobic'},6,1); repmat({'N-lim anaerobic'},6,1)];
names = unique(grp, 'stable');
fprintf('%-22s %7s %8s %8s | %7s %7s | %7s %7s', 'form', 'growth', 'glycerol', 'ethanol', 'fluxR2', 'fluxMRE', 'grR2', 'grMRE');
fprintf(' %12s', names{:}); fprintf('\n');
for k = 1:size(V,1)
    m = V{k,2};
    an = anaerobicModel(m);
    an = setParam(setParam(an, 'eq', 'r_1714', -23), 'obj', 'r_4041', 1);
    s = solveLP(an, 1);
    x = @(r) s.x(strcmp(an.rxns, r));
    figure('Visible','off'); [r2f, mre] = fluxPredSim(anaerobicModel(m)); close all
    figure('Visible','off'); [~, e, g] = growthNewSim(m); close all
    o = double(e(:,4)); p = g(:,4); rel = abs(p-o)./o;
    fprintf('%-22s %7.4f %8.3f %8.2f | %7.4f %7.4f | %7.4f %7.4f', V{k,1}, x('r_4041'), x('r_1808'), x('r_1761'), r2f, mre, ...
        1 - sum((o-p).^2)/sum((o-mean(o)).^2), mean(rel));
    for n = names', fprintf(' %12.4f', mean(rel(strcmp(grp, n{1})))); end
    fprintf('\n');
end
end

function gcvReverse()
% gcvReverse  Glycine cleavage r_0501 (glycine + NAD + THF => CH2-THF + NH4 + CO2 + NADH) with
%   the forward direction blocked (as released) and the reverse direction allowed
%   (CH2-THF + NH4 + CO2 + NADH => glycine + NAD + THF; lower bound -1000, upper bound 0),
%   and with both directions open. Maaheimo 2001 reports labelling from the glycine-forming
%   direction; the model's reactions are irreversible, so that direction cannot carry flux.
%   Sjoberg batch (glucose 23), 13C flux benchmark, growth benchmark.

paths = setupStudy('3.0.0-rc1');
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
base = loadTaggedModel(paths.yeast, '9.1.0');
rv = setParam(setParam(base, 'lb', 'r_0501', -1000), 'ub', 'r_0501', 0);
both = setParam(setParam(base, 'lb', 'r_0501', -1000), 'ub', {'r_0501','r_0507','r_0509'}, 1000);
V = {'forward blocked (released)', base; 'reverse only (r_0501)', rv; 'forward and reverse (r_0501)', both};
grp = [repmat({'N-lim aerobic'},9,1); repmat({'C-lim aerobic'},11,1); repmat({'C-lim anaerobic'},6,1); repmat({'N-lim anaerobic'},6,1)];
names = unique(grp, 'stable');
fprintf('%-30s %7s %8s %9s %9s | %7s %7s | %7s %7s', 'form', 'growth', 'glycerol', 'gly exp', 'r_0501', 'fluxR2', 'fluxMRE', 'grR2', 'grMRE');
fprintf(' %12s', names{:}); fprintf('\n');
for k = 1:size(V,1)
    m = V{k,2}; an = anaerobicModel(m);
    an = setParam(setParam(an, 'eq', 'r_1714', -23), 'obj', 'r_4041', 1);
    s = solveLP(an, 1); x = @(q) s.x(strcmp(an.rxns, q));
    figure('Visible','off'); [r2f, mre] = fluxPredSim(anaerobicModel(m)); close all
    figure('Visible','off'); [~, e, g] = growthNewSim(m); close all
    o = double(e(:,4)); p = g(:,4); rel = abs(p-o)./o;
    fprintf('%-30s %7.4f %8.3f %9.4f %9.4f | %7.4f %7.4f | %7.4f %7.4f', V{k,1}, x('r_4041'), x('r_1808'), x('r_1810'), x('r_0501'), r2f, mre, ...
        1 - sum((o-p).^2)/sum((o-mean(o)).^2), mean(rel));
    for n = names', fprintf(' %12.4f', mean(rel(strcmp(grp, n{1})))); end
    fprintf('\n');
end
end

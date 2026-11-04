function whatIfGrowthBench(traceCap)
% whatIfGrowthBench  Chemostat growth benchmark (Figure 3A-B, v9.1.0 test
%   code) for v9.1.0 with secretion of non-product exchanges capped at
%   traceCap, applied to the aerobic model before the benchmark adds its
%   own anaerobic and nitrogen-limitation changes.

paths = setupStudy('3.0.0-rc1');
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
base = loadTaggedModel(paths.yeast, '9.1.0');
[~, exch] = getExchangeRxns(base);
allowed = {'r_1761','r_1672','r_1808','r_1634','r_2056','r_2033','r_1546','r_1551', ...
           'r_2100','r_1832','r_4041','r_2111'};
capped = base;
closeIdx = exch(~ismember(base.rxns(exch), allowed));
capped.ub(closeIdx) = min(capped.ub(closeIdx), traceCap);
groups = [repmat({'N-limited aerobic'},9,1); repmat({'C-limited aerobic'},11,1);
          repmat({'C-limited anaerobic'},6,1); repmat({'N-limited anaerobic'},6,1)];
for mdl = {'unmodified', base; 'capped', capped}'
    figure('Visible','off'); [~, e, m] = growthNewSim(mdl{2}); close all
    o = double(e(:,4)); p = m(:,4);
    r2 = 1 - sum((o-p).^2)/sum((o-mean(o)).^2); rel = abs(p-o)./o;
    fprintf('%-10s R2 %.4f MRE %.4f |', mdl{1}, r2, mean(rel));
    for g = unique(groups)'
        fprintf(' %s %.4f', g{1}, mean(rel(strcmp(groups, g{1}))));
    end
    fprintf('\n');
end
end

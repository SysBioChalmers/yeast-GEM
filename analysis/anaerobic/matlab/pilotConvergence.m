function pilotConvergence(cap, thinnings, tag)
% pilotConvergence  How much thinning CHRR needs on v9.1.0 + new constraints,
%   Sjoberg batch (glucose 23, growth >= 99% of maximum): split-R-hat and
%   effective sample size for every reaction that is not fixed, 4 chains of
%   1000 samples, at two thinning levels.

paths = setupStudy('3.0.0-rc1');
out = fullfile(paths.results, 'pilot');
if ~isfolder(out), mkdir(out); end
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
model = anaerobicModel(loadTaggedModel(paths.yeast, '9.1.0'));
model = setParam(model, 'eq', 'r_1714', -23);
sol = solveLP(model);
model = setParam(model, 'lb', 'r_4041', 0.99 * sol.x(getIndexes(model, 'r_4041', 'rxns')));

watch = {'r_4041','r_1808','r_1761','r_1672','r_0713','r_0714','r_0658','r_0659', ...
         'r_0470','r_0471','r_0472','r_0226','r_0227','r_1654'};
res = table();
opts.cap = cap;
for thinning = thinnings
    [X, info] = chrrChains(model, 1000, thinning, 10*thinning, 4, 1, opts);
    fprintf('thinning %d: dimension %d, setup %.0f s, walking %.0f s, skipped steps %.3g%%\n', thinning, info.nDimensions, info.tSetup, info.tSteps, 100*max(info.skippedFraction));
    free = find(~ismember(model.rxns, info.fixedRxns));
    R = zeros(numel(free),1); E = R;
    for i = 1:numel(free)
        [R(i), E(i)] = mcmcDiag(squeeze(X(free(i),:,:)));
    end
    v = X(free,:,:); rngs = max(v,[],[2 3]) - min(v,[],[2 3]);
    fprintf('  sampled range: median %.3g, max %.3g; growth range %.4g\n', median(rngs), max(rngs), range(reshape(X(strcmp(model.rxns,'r_4041'),:,:),1,[])));
    fprintf('  all free reactions (%d): R-hat max %.3f, 99th pct %.3f; ESS min %.0f, 1st pct %.0f, median %.0f\n', ...
        numel(free), max(R), prctile(R,99), min(E), prctile(E,1), median(E));
    for w = watch
        k = find(strcmp(model.rxns, w{1}));
        [r, e] = mcmcDiag(squeeze(X(k,:,:)));
        v = X(k,:,:);
        res = [res; table(thinning, w, r, e, mean(v(:)), prctile(v(:),2.5), prctile(v(:),97.5), ...
            'VariableNames', {'thinning','rxn','rhat','ess','mean','q2_5','q97_5'})]; %#ok<AGROW>
    end
    save(fullfile(out, sprintf('conv_%s_thin%d.mat', tag, thinning)), 'X', 'info', '-v7.3');
end
writetable(res, fullfile(out, sprintf('convergence_watch_%s.tsv', tag)), 'FileType','text','Delimiter','\t');
disp(res);
end

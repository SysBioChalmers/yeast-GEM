function pilotRHMC(nSamples, cap)
% pilotRHMC  RHMC on v9.1.0 + new anaerobic constraints, Sjoberg batch
%   (glucose 23), growth >= 99% of maximum: run time, mixing, and whether the
%   growth range is explored; then the cycle-free projection of the samples.

paths = setupStudy('3.0.0-rc1');
addRHMC();
out = fullfile(paths.results, 'pilot');
if ~isfolder(out), mkdir(out); end
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
model = anaerobicModel(loadTaggedModel(paths.yeast, '9.1.0'));
model = setParam(model, 'eq', 'r_1714', -23);
sol = solveLP(model);
bio = getIndexes(model, 'r_4041', 'rxns');
gmax = sol.x(bio);
model = setParam(model, 'lb', 'r_4041', 0.99 * gmax);

[X, dg] = rhmcSample(model, nSamples, struct('cap', cap, 'seed', 1));
u = dg.Properties.UserData;
moving = dg.max - dg.min > 1e-9;
fprintf('RHMC: %d draws in %.0f s; %d reactions vary\n', u.draws, u.seconds, sum(moving));
fprintf('  R-hat max %.3f, 99th pct %.3f; ESS min %.0f, 1st pct %.0f, median %.0f (varying reactions)\n', ...
    max(dg.rhat(moving)), prctile(dg.rhat(moving),99), min(dg.ess(moving)), prctile(dg.ess(moving),1), median(dg.ess(moving)));
fprintf('  growth sampled %.4f-%.4f (feasible %.4f-%.4f)\n', min(X(bio,:)), max(X(bio,:)), 0.99*gmax, gmax);

[~, exch] = getExchangeRxns(model);
tic; W = cycleFree(model, X, [exch(:); bio]); tCF = toc;
internal = true(numel(model.rxns),1); internal(exch) = false;
fprintf('cycle-free projection: %.0f s; max |internal flux| before %.1f, after %.1f (cap %g)\n', ...
    tCF, max(abs(X(internal,:)),[],'all'), max(abs(W(internal,:)),[],'all'), cap);
save(fullfile(out, sprintf('rhmc_cap%g_n%d.mat', cap, nSamples)), 'X', 'W', 'dg', '-v7.3');
writetable(dg, fullfile(out, sprintf('rhmc_cap%g_n%d_diag.tsv', cap, nSamples)), 'FileType','text','Delimiter','\t');
end

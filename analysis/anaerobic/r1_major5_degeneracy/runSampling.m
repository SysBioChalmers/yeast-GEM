function runSampling(variant, scenario, nSamples, cap, seed, optFrac)
% runSampling  Sample one (variant, scenario) with RHMC, remove cycles, save.
%   runSampling('v910_new', 'sjoberg', 1000, 200, 1)
%   runSampling('v910_new', 'sjoberg', 1000, 200, 1, 0.99)   % sensitivity case
%
%   Writes results/sampling/<variant>__<scenario>.mat with the raw samples X,
%   the cycle-free samples W, the per-reaction mixing diagnostics, the pFBA
%   solution and the scenario spec. Each job writes only its own file.

paths = setupStudy('3.0.0-rc1');
addRHMC();
out = fullfile(paths.results, 'sampling');
if ~isfolder(out), mkdir(out); end
% anaerobicModel (original constraints) reads its data relative to code/otherChanges.
cd(fullfile(paths.yeast, 'code', 'otherChanges'));

if nargin < 6, optFrac = 0.99999; end
[model, spec] = scenarioModel(paths, variant, scenario, optFrac);
file = sprintf('%s__%s__opt%g', variant, scenario, optFrac);
if seed ~= 1, file = sprintf('%s__seed%d', file, seed); end   % first chain keeps the original name
if cap ~= 200, file = sprintf('%s__cap%g', file, cap); end
file = fullfile(out, [file '.mat']);
rxns = model.rxns;
if ~spec.feasible
    fprintf('%s / %s: infeasible (no growth); nothing to sample\n', variant, scenario);
    save(file, 'spec', 'rxns');
    return
end

[X, dg] = rhmcSample(model, nSamples, struct('cap', cap, 'seed', seed));
[~, exch] = getExchangeRxns(model);
bio = getIndexes(model, 'r_4041', 'rxns');
W = cycleFree(model, X, [exch(:); bio]);

internal = true(numel(rxns),1); internal(exch) = false;
spec.maxInternalCycleFree = max(abs(W(internal,:)), [], 'all');
spec.cap = cap; spec.seed = seed;
moving = dg.max - dg.min > 1e-9;
u = dg.Properties.UserData;
fprintf('%s / %s: %d draws in %.0f s; R-hat max %.3f; ESS min %.0f median %.0f; max cycle-free internal flux %.1f (cap %g)\n', ...
    variant, scenario, u.draws, u.seconds, max(dg.rhat(moving)), min(dg.ess(moving)), median(dg.ess(moving)), ...
    spec.maxInternalCycleFree, cap);
save(file, 'X', 'W', 'dg', 'spec', 'rxns', '-v7.3');
end

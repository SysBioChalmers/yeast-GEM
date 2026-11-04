function pilotSampling()
% pilotSampling  Size the sampling problem on v9.1.0 + new anaerobic
%   constraints, Sjoberg batch (glucose 23): loop reactions, polytope
%   dimension, and time per sample.

paths = setupStudy('3.0.0-rc1');
out = fullfile(paths.results, 'pilot');
if ~isfolder(out), mkdir(out); end
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
model = anaerobicModel(loadTaggedModel(paths.yeast, '9.1.0'));
model = setParam(model, 'eq', 'r_1714', -23);

% Reactions that can carry flux with every exchange closed are in internal
% cycles (loops).
[~, exch] = getExchangeRxns(model);
closed = model;
closed.lb(exch) = 0; closed.ub(exch) = 0;
% Release forced fluxes (e.g. the non-growth maintenance reaction), otherwise
% the closed model is infeasible and every range comes back empty.
fprintf('reactions with forced flux: %s\n', strjoin(model.rxns(model.lb > 0 | model.ub < 0)', ' '));
closed.lb(closed.lb > 0) = 0; closed.ub(closed.ub < 0) = 0;
assert(~isempty(solveLP(closed).x), 'closed model infeasible');
tic; [lo, hi] = getAllowedBounds(closed); tLoop = toc;
loopRxns = model.rxns(abs(lo) > 1e-6 | abs(hi) > 1e-6);
fprintf('exchange reactions: %d; loop reactions: %d (FVA %.0f s)\n', numel(exch), numel(loopRxns), tLoop);
writetable(table(loopRxns(:), constructEquations(model, loopRxns(:)), ...
    'VariableNames', {'rxn','equation'}), fullfile(out, 'loop_rxns_v910_new.tsv'), ...
    'FileType','text','Delimiter','\t');

% Near-optimal space: growth >= 99% of the maximum.
sol = solveLP(model);
bio = getIndexes(model, 'r_4041', 'rxns');
fprintf('max growth %.4f\n', -sol.f);
m99 = setParam(model, 'lb', 'r_4041', 0.99 * sol.x(bio));

tic; [s, info] = sampleCHRR(m99, 50, 100, 100, 1); tSample = toc;
fprintf('dimension %d, fixed reactions %d, MVE converged %d\n', info.nDimensions, numel(info.fixedRxns), info.mveConverged);
fprintf('50 samples x 100 thinning + 100 burn-in: %.0f s total\n', tSample);
save(fullfile(out, 'pilot.mat'), 's', 'info', 'loopRxns', 'tLoop', 'tSample');
end

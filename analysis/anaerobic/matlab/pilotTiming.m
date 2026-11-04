function pilotTiming()
% pilotTiming  Separate sampleCHRR's one-off setup cost (FVA, rounding) from
%   its cost per hit-and-run step, on v9.1.0 + new constraints, Sjoberg batch.

paths = setupStudy('3.0.0-rc1');
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
model = anaerobicModel(loadTaggedModel(paths.yeast, '9.1.0'));
model = setParam(model, 'eq', 'r_1714', -23);
sol = solveLP(model);
m99 = setParam(model, 'lb', 'r_4041', 0.99 * sol.x(getIndexes(model, 'r_4041', 'rxns')));

tic; sampleCHRR(m99, 10, 10, 0, 1); t1 = toc;          % 100 steps
tic; sampleCHRR(m99, 10, 20000, 0, 1); t2 = toc;       % 200,000 steps
perStep = (t2 - t1) / (200000 - 100);
fprintf('setup %.0f s, per step %.2e s, so 1e6 steps take %.0f s\n', t1, perStep, 1e6*perStep);
end

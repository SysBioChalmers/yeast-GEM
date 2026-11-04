function pilotWorst()
% pilotWorst  Which reactions mix worst in the convergence pilot, and are
%   they the loop reactions?

paths = setupStudy('3.0.0-rc1');
out = fullfile(paths.results, 'pilot');
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
model = anaerobicModel(loadTaggedModel(paths.yeast, '9.1.0'));
loops = readtable(fullfile(out, 'loop_rxns_v910_new.tsv'), 'FileType','text','Delimiter','\t');
isLoop = ismember(model.rxns, loops.rxn);
load(fullfile(out, 'conv_thin10000.mat'), 'X', 'info');
free = find(~ismember(model.rxns, info.fixedRxns));
R = zeros(numel(free),1); E = R; rng_ = R;
for i = 1:numel(free)
    [R(i), E(i)] = mcmcDiag(squeeze(X(free(i),:,:)));
    v = X(free(i),:,:); rng_(i) = max(v(:)) - min(v(:));
end
fprintf('free reactions: %d, of which loop reactions: %d\n', numel(free), sum(isLoop(free)));
for grp = {'loop', 'not loop'}
    sel = isLoop(free) == strcmp(grp{1}, 'loop');
    fprintf('%-9s R-hat median %.3f max %.3f | ESS median %.0f min %.0f | sampled range median %.1f max %.1f\n', ...
        grp{1}, median(R(sel)), max(R(sel)), median(E(sel)), min(E(sel)), median(rng_(sel)), max(rng_(sel)));
end
[~, o] = sort(R, 'descend');
T = table(model.rxns(free(o(1:25))), model.rxnNames(free(o(1:25))), isLoop(free(o(1:25))), R(o(1:25)), E(o(1:25)), rng_(o(1:25)), ...
    'VariableNames', {'rxn','name','loop','rhat','ess','range'});
disp(T);
writetable(T, fullfile(out, 'worst_mixing.tsv'), 'FileType','text','Delimiter','\t');
end

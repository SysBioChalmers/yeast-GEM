function fvaRun(variant, optFrac)
% fvaRun  Flux variability analysis for one model variant, all 13 scenarios of
%   the sampling (scenarioModel), with the same constraints (optFrac as there:
%   0.99999 optimal face, 1.01 a +-1% band on glucose uptake and growth).
%   Reactions: measured exchanges, the reactions discussed in the text, and the
%   reactions compared in the fluxomics benchmark. Model bounds are used as
%   they are (no internal cap). Each extreme is also projected onto a loop-free
%   flux vector with the exchanges and growth kept (cycleFree), which gives the
%   loop-free value of that reaction (CycleFreeFlux, Desouki et al. 2015); the
%   raw extreme of a reaction in a loop is otherwise set by its +-1000 bound.
%   Writes results/fva/<variant>__opt<optFrac>.tsv.

paths = setupStudy('3.0.0-rc1');
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
out = fullfile(paths.results, 'fva');
if ~isfolder(out), mkdir(out); end
F = readtable(fullfile(paths.yeast, 'data', 'physiology', 'flux_data_anaerobic.tsv'), ...
    'FileType','text','ReadVariableNames',false,'Delimiter','\t');
base = {'r_4041','r_1714','r_1808','r_1761','r_1672','r_1634','r_2056','r_2033', ...
        'r_0713','r_0714','r_2131','r_0659','r_0658','r_0471','r_0470','r_0476','r_0472', ...
        'r_0501','r_0487','r_0226','r_0227','r_1654','r_1115','r_0491','r_0252','r_1258','r_4224', ...
        'r_1824','r_2129','r_1810'};
scen = {'sjoberg','bakker_wt','bakker_adh3','bakker_adh3adh4','boles_wt','boles_mae1'};
sets = unique(F{:,6});
for s = sets', scen{end+1} = ['flux_' s{1}]; end %#ok<AGROW>

T = table();
for sc = scen
    [m, spec] = scenarioModel(paths, variant, sc{1}, optFrac);
    if ~spec.feasible, continue; end
    rx = base;
    if startsWith(sc{1}, 'flux_')
        rx = union(rx, F{strcmp(F{:,6}, extractAfter(sc{1}, 'flux_')), 7}, 'stable');
    end
    rx = rx(ismember(rx, m.rxns) & ~cellfun(@isempty, rx));
    [~, exch] = getExchangeRxns(m);
    keep = [exch(:); getIndexes(m, 'r_4041', 'rxns')];
    for k = 1:numel(rx)
        i = find(strcmp(m.rxns, rx{k}));
        v = NaN(1,4);
        for dirn = [-1 1]
            s = solveLP(setParam(m, 'obj', rx{k}, dirn));
            if isempty(s.x), continue; end
            w = cycleFree(m, s.x, keep);
            col = (dirn + 3) / 2;               % 1 = min, 2 = max
            v(col) = s.x(i); v(col + 2) = w(i);
        end
        pf = NaN; if ~isempty(spec.pfba), pf = spec.pfba(strcmp(m.rxns, rx{k})); if isempty(pf), pf = NaN; end, end
        T = [T; table({variant}, sc, rx(k), m.rxnNames(i), pf, v(1), v(2), v(3), v(4), ...
            'VariableNames', {'variant','scenario','rxn','name','pfba','min','max','min_loopfree','max_loopfree'})]; %#ok<AGROW>
    end
    fprintf('%s / %s: %d reactions\n', variant, sc{1}, numel(rx));
end
writetable(T, fullfile(out, sprintf('%s__opt%g.tsv', variant, optFrac)), 'FileType','text','Delimiter','\t');
end

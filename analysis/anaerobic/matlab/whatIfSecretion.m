function whatIfSecretion(traceCap)
% whatIfSecretion  What-if analysis, no change to the released model: close
%   secretion of every exchange that is not a known anaerobic product to
%   traceCap mmol/gDW/h (the wild type needs trace byproducts such as glycine
%   and glycolaldehyde to leave the cell), and
%   report (a) whether gpd1 gpd2 can still grow anaerobically and through
%   which exit, (b) the effect on the wild-type predictions (Sjoberg,
%   Bakker wild type and adh3) and the fluxomics benchmark.

paths = setupStudy('3.0.0-rc1');
out = fullfile(paths.results, 'whatif_secretion');
if ~isfolder(out), mkdir(out); end
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
base = anaerobicModel(loadTaggedModel(paths.yeast, '9.1.0'));
[~, exch] = getExchangeRxns(base);

% What does the unmodified model secrete (pFBA, Sjoberg and Bakker batch)?
secreted = false(numel(base.rxns),1);
for g = [23 19.1]
    s = solveLP(setParam(base, 'eq', 'r_1714', -g), 1);
    secreted(exch(s.x(exch) > 1e-6)) = true;
end
fprintf('Exchanges secreted by the unmodified v9.1.0/new model (pFBA, glucose 23 and 19.1):\n');
for k = find(secreted)'
    fprintf('  %-8s %s\n', base.rxns{k}, base.rxnNames{k});
end

% Known products of anaerobic S. cerevisiae on glucose, plus the inorganic
% outputs. Everything else may still be taken up but not secreted.
allowed = {'r_1761','r_1672','r_1808','r_1634','r_2056','r_2033','r_1546','r_1551', ...  % ethanol CO2 glycerol acetate succinate pyruvate lactate(D/L)
           'r_2100','r_1832','r_4041','r_2111'};                                         % water H+ biomass growth
[inModel, pos] = ismember(allowed, base.rxns);
fprintf('\nAllowed secretions: %s\n', strjoin(strcat(allowed(inModel), ' (', base.rxnNames(pos(inModel))', ')'), ', '));
fprintf('Other secretions capped at %g\n', traceCap);
capped = base;
closeIdx = exch(~ismember(base.rxns(exch), allowed));
capped.ub(closeIdx) = min(capped.ub(closeIdx), traceCap);

T = table();
for mdl = {'unmodified', base; sprintf('secretion capped %g', traceCap), capped}'
    m = mdl{2};
    for sc = {'WT', {}, 19.1; 'gpd1 gpd2', {'YDL022W','YOL059W'}, 19.1; 'adh3', {'YMR083W'}, 12.6; 'WT (Sjoberg)', {}, 23}'
        mm = m;
        if ~isempty(sc{2}), mm = removeGenes(mm, sc{2}); end
        mm = setParam(mm, 'eq', 'r_1714', -sc{3});
        mm = setParam(mm, 'obj', 'r_4041', 1);
        s = solveLP(mm, 1);
        if isempty(s.x), mu = 0; gly = NaN; eth = NaN; else
            mu = s.x(strcmp(mm.rxns,'r_4041')); gly = s.x(strcmp(mm.rxns,'r_1808')); eth = s.x(strcmp(mm.rxns,'r_1761')); end
        T = [T; table(mdl(1), sc(1), sc{3}, mu, gly, eth, 'VariableNames', {'model','strain','glucose','growth','glycerol','ethanol'})]; %#ok<AGROW>
        if strcmp(sc{1}, 'gpd1 gpd2') && ~isempty(s.x) && mu > 1e-6
            fprintf('\n%s, gpd1 gpd2 still grows (%.4f); secreted exchanges:\n', mdl{1}, mu);
            for k = exch(s.x(exch) > 1e-3)'
                fprintf('  %-8s %-35s %.3f\n', mm.rxns{k}, mm.rxnNames{k}, s.x(k));
            end
        end
    end
end
disp(T);
writetable(T, fullfile(out, sprintf('growth_cap%g.tsv', traceCap)), 'FileType','text','Delimiter','\t');

% Fluxomics benchmark with the secretion closed.
figure('Visible','off');
[R2a, mreA] = fluxPredSim(base); close all
figure('Visible','off');
[R2b, mreB] = fluxPredSim(capped); close all
fprintf('\nFluxomics benchmark: unmodified R2 %.4f MRE %.4f | secretion capped R2 %.4f MRE %.4f\n', R2a, mreA, R2b, mreB);
end

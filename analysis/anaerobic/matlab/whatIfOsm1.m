function whatIfOsm1()
% whatIfOsm1  What-if (no change to the released model): fumarate reductase
%   curation from the literature review.
%     r_0455 (cytosolic FADH2 + fumarate) gene rule YEL047C or YJR051W, irreversible
%     r_0454 (mitochondrial FADH2 + fumarate) gene rule YJR051W, irreversible
%     r_1000 (cytosolic FMNH2 + fumarate) gene rule YEL047C or YJR051W, irreversible
%     r_4264 (mitochondrial NADH + fumarate, OSM1) removed (bounds 0)
%   Scored on frd1, osm1, frd1 osm1 and gpd1 gpd2 anaerobic growth, the
%   wild-type Sjoberg and Bakker predictions, the fluxomics and growth benchmarks.

paths = setupStudy('3.0.0-rc1');
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
for ver = {'9.1.0', '9.1.1'}
    base = loadTaggedModel(paths.yeast, ver{1});
    cur = changeGrRules(base, 'r_0455', 'YEL047C or YJR051W', true);
    cur = changeGrRules(cur, 'r_0454', 'YJR051W', true);
    cur = changeGrRules(cur, 'r_1000', 'YEL047C or YJR051W', true);
    cur = setParam(cur, 'lb', {'r_0455','r_0454','r_1000'}, 0);
    cur = setParam(cur, 'eq', 'r_4264', 0);
    for mdl = {'released', base; 'Osm1 curation', cur}'
        an = anaerobicModel(mdl{2});
        s0 = solveLP(setParam(setParam(an, 'eq', 'r_1714', -23), 'obj', 'r_4041', 1), 1);
        fprintf('\n%s %s: r_4264 flux (Sjoberg pFBA) %.4f\n', ver{1}, mdl{1}, s0.x(strcmp(an.rxns,'r_4264')));
        for sc = {'Sjoberg', {}, 23; 'Bakker WT', {}, 19.1; 'Bakker adh3', {'YMR083W'}, 12.6; ...
                  'frd1', {'YEL047C'}, 23; 'osm1', {'YJR051W'}, 23; 'frd1 osm1', {'YEL047C','YJR051W'}, 23; ...
                  'gpd1 gpd2', {'YDL022W','YOL059W'}, 19.1}'
            m = an; if ~isempty(sc{2}), m = removeGenes(m, sc{2}); end
            m = setParam(setParam(m, 'eq', 'r_1714', -sc{3}), 'obj', 'r_4041', 1);
            s = solveLP(m, 1);
            if isempty(s.x), fprintf('  %-12s infeasible\n', sc{1}); continue; end
            fprintf('  %-12s growth %.4f glycerol %.3f succinate %.3f\n', sc{1}, s.x(strcmp(m.rxns,'r_4041')), ...
                s.x(strcmp(m.rxns,'r_1808')), s.x(strcmp(m.rxns,'r_2056')));
        end
        figure('Visible','off'); [R2f, mref] = fluxPredSim(an); close all
        figure('Visible','off'); [R2g, e, g] = growthNewSim(mdl{2}); close all
        fprintf('  fluxomics R2 %.4f MRE %.4f; growth R2 %.4f MRE %.4f\n', R2f, mref, R2g, ...
            mean(abs(g(:,4)-double(e(:,4)))./double(e(:,4))));
    end
end
end

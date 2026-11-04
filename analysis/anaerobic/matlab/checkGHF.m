function checkGHF()
% checkGHF  Checks of package G+H+F (glycine cleavage open; Hom6 NADPH only;
%   FMNH2-linked fumarate reductase r_1000 removed) against v9.1.0:
%   1. gpd1 gpd2 anaerobic growth, and its rescue by acetaldehyde or acetoin
%      supplied in the medium (Ansell 1997; Bjorkqvist 1997).
%   2. Anaerobic single-gene deletions (Sjoberg batch): genes whose predicted
%      growth ratio differs between the released model and the package.
%   3. Aerobic gene essentiality (yeast-GEM essentialGenes test).

paths = setupStudy('3.0.0-rc1');
out = fullfile(paths.results, 'whatif_c1');
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
base = loadTaggedModel(paths.yeast, '9.1.0');
pk = setParam(base, 'ub', {'r_0501','r_0507','r_0509','r_0546'}, [1000 1000 1000 0]);
pk = setParam(setParam(pk, 'lb', 'r_1000', 0), 'ub', 'r_1000', 0);
models = {'released', base; 'G+H+F', pk};

acetoinE = find(strcmp(base.metNames, '(R)-acetoin') & strcmp(base.comps(base.metComps), 'e'));
fprintf('(R)-acetoin[e] present: %d\n', ~isempty(acetoinE));
for k = 1:2
    an = anaerobicModel(models{k,2});
    if ~isempty(acetoinE)
        r.rxns = {'acetoin_uptake'}; r.mets = an.mets(acetoinE); r.stoichCoeffs = 1; r.lb = 0; r.ub = 0;
        r.rxnNames = {'(R)-acetoin uptake (test)'};
        an = addRxns(an, r, 1, [], false);
    end
    ko = removeGenes(an, {'YDL022W','YOL059W'});
    for sup = {'none', {}, 0; 'acetaldehyde 2', {'r_1631'}, -2; 'acetaldehyde 5', {'r_1631'}, -5; ...
               'acetoin 2', {'acetoin_uptake'}, 2; 'acetoin 5', {'acetoin_uptake'}, 5}'
        for strain = {'WT', an; 'gpd1 gpd2', ko}'
            mm = setParam(setParam(strain{2}, 'eq', 'r_1714', -19.1), 'obj', 'r_4041', 1);
            if ~isempty(sup{2})
                if sup{3} < 0, mm = setParam(mm, 'lb', sup{2}, sup{3}); else, mm = setParam(mm, 'ub', sup{2}, sup{3}); end
            end
            s = solveLP(mm);
            g = 0; if ~isempty(s.x), g = s.x(strcmp(mm.rxns,'r_4041')); end
            fprintf('%-9s %-10s supplement %-15s growth %.4f\n', models{k,1}, strain{1}, sup{1}, g);
        end
    end
end

fprintf('\nAnaerobic single-gene deletions (glucose 23):\n');
gr = cell(1,2);
for k = 1:2
    an = setParam(setParam(anaerobicModel(models{k,2}), 'eq', 'r_1714', -23), 'obj', 'r_4041', 1);
    [genes, ~, ~, ~, ratio] = findGeneDeletions(an, 'sgd', 'fba');
    v = ones(numel(an.genes),1); v(genes) = ratio; gr{k} = containers.Map(an.genes, num2cell(v));
end
g = keys(gr{1}); n = 0;
for i = 1:numel(g)
    a = gr{1}(g{i}); b = gr{2}(g{i});
    if abs(a - b) > 0.02
        fprintf('  %-10s released %.3f  G+H+F %.3f\n', g{i}, a, b); n = n + 1;
    end
end
fprintf('  %d genes differ by more than 0.02 in growth ratio\n', n);

cd(fullfile(paths.yeast, 'code', 'modelTests'));
for k = 1:2
    [acc, tp, tn, fn, fp] = essentialGenes(models{k,2});
    fprintf('Aerobic essentiality %-9s accuracy %.4f (tp %d tn %d fn %d fp %d)\n', models{k,1}, acc, tp, tn, fn, fp);
end
end

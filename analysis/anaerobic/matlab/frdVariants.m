function frdVariants()
% frdVariants  The pieces of the fumarate reductase curation tested one at a time and together
%   (v9.1.0, updated anaerobic constraints): gene rules, irreversibility, OSM1 on r_0454, removal of
%   r_4264 (bounds 0 as the equivalent). Anaerobic batch, glucose 23: growth of wild type, frd1, osm1,
%   frd1 osm1 and gpd1 gpd2, Sjoberg glycerol. Then the flux range of the four reactions over all
%   optimal solutions of the released model (growth fixed at 99.999% of the optimum).

paths = setupStudy('3.0.0-rc1');
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
base = loadTaggedModel(paths.yeast, '9.1.0');
rules = @(m) changeGrRules(changeGrRules(m, 'r_0455', 'YEL047C or YJR051W', true), 'r_1000', 'YEL047C or YJR051W', true);
irrev = @(m) setParam(m, 'lb', {'r_0454','r_0455','r_1000'}, 0);
r0454 = @(m) changeGrRules(m, 'r_0454', 'YJR051W', true);
rm4264 = @(m) setParam(m, 'eq', 'r_4264', 0);
V = {'released', base; 'rules r_0455, r_1000 only', rules(base); 'irreversible only', irrev(base); 'r_4264 removed only', rm4264(base); ...
     'r_0454 OSM1 only', r0454(base); 'rules + irreversible', irrev(rules(base)); ...
     'rules + irreversible + r_0454', r0454(irrev(rules(base))); 'all (branch)', rm4264(r0454(irrev(rules(base))))};
ko = {'wild type', {}; 'frd1', {'YEL047C'}; 'osm1', {'YJR051W'}; 'frd1 osm1', {'YEL047C','YJR051W'}; 'gpd1 gpd2', {'YDL022W','YOL059W'}};
fprintf('%-34s', 'variant'); fprintf(' %10s', ko{:,1}); fprintf(' %9s\n', 'glycerol');
for k = 1:size(V,1)
    an = anaerobicModel(V{k,2});
    fprintf('%-34s', V{k,1});
    for j = 1:size(ko,1)
        m = an; if ~isempty(ko{j,2}), m = removeGenes(m, ko{j,2}); end
        m = setParam(setParam(m, 'eq', 'r_1714', -23), 'obj', 'r_4041', 1);
        s = solveLP(m, 1); g = 0; if ~isempty(s.x), g = s.x(strcmp(m.rxns, 'r_4041')); end
        fprintf(' %10.4f', g);
        if j == 1, gl = s.x(strcmp(m.rxns, 'r_1808')); end
    end
    fprintf(' %9.3f\n', gl);
end
an = setParam(setParam(anaerobicModel(base), 'eq', 'r_1714', -23), 'obj', 'r_4041', 1);
s = solveLP(an, 1); mu = s.x(strcmp(an.rxns, 'r_4041'));
fprintf('\nFlux ranges, released model, Sjoberg, growth fixed at 99.999%% of optimum (pFBA value, min, max):\n');
mm = setParam(an, 'lb', 'r_4041', 0.99999 * mu);
for r = {'r_0454','r_0455','r_1000','r_4264'}
    lo = solveLP(setParam(mm, 'obj', r{1}, -1)); hi = solveLP(setParam(mm, 'obj', r{1}, 1)); i = strcmp(an.rxns, r{1});
    fprintf('  %s  pFBA %8.4f  min %8.4f  max %8.4f\n', r{1}, s.x(i), lo.x(i), hi.x(i));
end
end

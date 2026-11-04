function geneDeletionScreen(variant)
% geneDeletionScreen  Single-gene deletion growth ratios for every gene of one
%   model variant (see variantModel), in four conditions:
%     minimal glucose medium, aerobic and anaerobic (anaerobic function of the
%       variant applied to the aerobic model);
%     synthetic complete medium (completeY7, as the yeast-GEM essentiality
%       test), aerobic and anaerobic.
%   Glucose uptake up to 10 (minimal) or 20 (complete) mmol/gDW/h, growth
%   maximised. Writes results/ko_screen/<variant>.tsv: gene and growth ratio
%   (mutant/wild type) per condition, plus the wild-type growth rates.

paths = setupStudy('3.0.0-rc1');
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
out = fullfile(paths.results, 'ko_screen');
if ~isfolder(out), mkdir(out); end
[~, base, anaerobicFn, post] = variantModel(paths, variant);
minimal = setParam(base, 'lb', 'r_1714', -10);
rich = completeY7(base);
conds = {'aer_minimal', minimal; 'an_minimal', post(anaerobicFn(minimal)); ...
         'aer_complete', rich; 'an_complete', post(anaerobicFn(rich))};
T = table(base.genes, 'VariableNames', {'gene'});
wt = zeros(1, size(conds,1));
for c = 1:size(conds,1)
    m = setParam(conds{c,2}, 'obj', 'r_4041', 1);
    s = solveLP(m); wt(c) = s.f * -1;
    if isempty(s.x), wt(c) = 0; else, wt(c) = s.x(strcmp(m.rxns, 'r_4041')); end
    [genes, ~, ~, ~, ratio] = findGeneDeletions(m, 'sgd', 'fba');
    r = ones(numel(m.genes), 1); r(genes) = ratio;
    [~, pos] = ismember(T.gene, m.genes);
    v = NaN(height(T), 1); v(pos > 0) = r(pos(pos > 0));
    T.(conds{c,1}) = v;
    fprintf('%s %s: wild-type growth %.4f; %d genes lethal (ratio < 1e-6)\n', variant, conds{c,1}, wt(c), nnz(v < 1e-6));
end
T.Properties.UserData = wt;
writetable(T, fullfile(out, [variant '.tsv']), 'FileType','text','Delimiter','\t');
writetable(cell2table([conds(:,1)'; num2cell(wt)]), fullfile(out, [variant '_wildtype.tsv']), 'FileType','text','Delimiter','\t', 'WriteVariableNames', false);
end

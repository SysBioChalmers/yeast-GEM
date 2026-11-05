function gpdViability(variants, outName)
% gpdViability  Predicted anaerobic growth of glycerol-3-phosphate
%   dehydrogenase mutants (gpd1, gpd2, gpd1 gpd2) in each model variant,
%   anaerobic batch with glucose uptake fixed at 19.1 (Bakker 2000 wild type).
%   Literature: gpd1 gpd2 cannot grow anaerobically and gpd2 grows poorly
%   (Ansell et al. 1997 EMBO J; Bjorkqvist et al. 1997 AEM) -- to be checked
%   against the full texts.

paths = setupStudy('3.0.0-rc1');
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
if nargin < 1, variants = {'vdev_new','vdev_old','v902_new','v902_old'}; end
if nargin < 2, outName = 'gpd_viability.tsv'; end
out = fullfile(paths.results, outName);
kos = {'none', {}; 'gpd1', {'YDL022W'}; 'gpd2', {'YOL059W'}; 'gpd1 gpd2', {'YDL022W','YOL059W'}};
T = table();
for v = variants
    m0 = variantModel(paths, v{1});
    for k = 1:size(kos,1)
        m = m0;
        if ~isempty(kos{k,2}), m = removeGenes(m, kos{k,2}); end
        m = setParam(m, 'eq', 'r_1714', -19.1);
        m = setParam(m, 'obj', 'r_4041', 1);
        sol = solveLP(m, 1);
        if isempty(sol.x)
            mu = 0; gly = NaN;
        else
            mu = sol.x(strcmp(m.rxns, 'r_4041')); gly = sol.x(strcmp(m.rxns, 'r_1808'));
        end
        T = [T; table(v, kos(k,1), mu, gly, 'VariableNames', {'variant','knockout','growth','glycerol'})]; %#ok<AGROW>
    end
end
writetable(T, out, 'FileType','text','Delimiter','\t');
disp(T);
end

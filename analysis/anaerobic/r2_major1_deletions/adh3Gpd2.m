function adh3Gpd2()
% adh3Gpd2  Anaerobic batch predictions (pFBA) for adh3, gpd2 and adh3 gpd2,
%   as in knockoutPanel (glucose uptake fixed at 20 mmol/gDW/h, growth
%   maximised). Writes results/adh3_gpd2.tsv.
variants = {'vdev_new','vdev_old','v902_new','v902_old'};
paths = setupStudy('3.0.0-rc1');
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
strains = {'wild type',{}; 'adh3',{'YMR083W'}; 'gpd2',{'YOL059W'}; 'adh3 gpd2',{'YMR083W','YOL059W'}; 'mae1',{'YKL029C'}};
glc = 20; mw = 180.156; T = table();
for v = variants
    an = variantModel(paths, v{1});
    for k = 1:size(strains,1)
        ko = strains{k,2};
        m = an; if ~isempty(ko), m = removeGenes(m, ko(ismember(ko, an.genes))); end
        m = setParam(setParam(m,'eq','r_1714',-glc),'obj','r_4041',1);
        s = solveLP(m,1);
        if isempty(s.x), mu = 0; else, mu = s.x(strcmp(m.rxns,'r_4041')); end
        if mu < 1e-6, vals = [0 0 NaN NaN]; else
            x = @(r) s.x(strcmp(m.rxns,r));
            vals = [mu, mu/(glc*mw/1000), x('r_1761')*46.07/(glc*mw), x('r_1808')*92.09/(glc*mw)];
        end
        T = [T; table(string(v{1}), string(strains{k,1}), vals(1), vals(2), vals(3), vals(4), ...
            'VariableNames', {'variant','strain','mu','Yx','Yethanol','Yglycerol'})]; %#ok<AGROW>
    end
end
writetable(T, fullfile(paths.results, 'adh3_gpd2.tsv'), 'FileType','text','Delimiter','\t');
disp(T)

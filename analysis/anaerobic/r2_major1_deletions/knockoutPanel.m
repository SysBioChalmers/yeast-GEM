function knockoutPanel(variants, outName)
% knockoutPanel  Anaerobic batch predictions (pFBA) for every deletion strain
%   with anaerobic physiology data (data/knockout_data.tsv and
%   data/knockout_data_extended.tsv) and the fumarate reductase mutants.
%   Minimal glucose medium with the variant's anaerobic function; glucose
%   uptake fixed at 20 mmol/gDW/h (yields are compared, which barely depend on
%   the uptake rate); growth maximised. Writes results/ko_panel.tsv with
%   growth, biomass/ethanol/glycerol yields (g/g glucose), acetate and
%   succinate (mmol/g biomass) per variant and strain.

if nargin < 1, variants = {'vdev_new','vdev_old','v902_new','v902_old'}; end
if nargin < 2, outName = 'ko_panel.tsv'; end
paths = setupStudy('3.0.0-rc1');
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
strains = {
    'wild type',        {}
    'adh3',             {'YMR083W'}
    'mae1',             {'YKL029C'}
    'ald6',             {'YPL061W'}
    'ald5',             {'YER073W'}
    'ald4',             {'YOR374W'}
    'ald6 ald4 ald5',   {'YPL061W','YOR374W','YER073W'}
    'kgd1',             {'YIL125W'}
    'gdh1',             {'YOR375C'}
    'gdh1 gdh3',        {'YOR375C','YAL062W'}
    'sdh1',             {'YKL148C'}
    'gpd1',             {'YDL022W'}
    'gpd2',             {'YOL059W'}
    'gpd1 gpd2',        {'YDL022W','YOL059W'}
    'gpp1',             {'YIL053W'}
    'gpp2',             {'YER062C'}
    'gpp1 gpp2',        {'YIL053W','YER062C'}
    'frd1',             {'YEL047C'}
    'osm1',             {'YJR051W'}
    'frd1 osm1',        {'YEL047C','YJR051W'}
    'ole1',             {'YGL055W'}
    'plb2',             {'YMR006C'}
    'izh2',             {'YOL002C'}
    'ina1 (ylr413w)',   {'YLR413W'}
    'eug1',             {'YDR518W'}
    'yor012w',          {'YOR012W','YOR013W'}
    };
glc = 20; mw = 180.156;
T = table();
for v = variants
    an = variantModel(paths, v{1});
    for k = 1:size(strains,1)
        ko = strains{k,2};
        inModel = ismember(ko, an.genes);
        m = an; if any(inModel), m = removeGenes(m, ko(inModel)); end
        m = setParam(setParam(m, 'eq', 'r_1714', -glc), 'obj', 'r_4041', 1);
        s = solveLP(m, 1);
        x = @(r) s.x(strcmp(m.rxns, r));
        if isempty(s.x) || x('r_4041') < 1e-6
            vals = [0 0 NaN NaN NaN NaN NaN NaN];
        else
            mu = x('r_4041');
            vals = [mu, mu/(glc*mw/1000), x('r_1761')*46.07/(glc*mw), x('r_1808')*92.09/(glc*mw), ...
                    x('r_1634')/mu, x('r_2056')/mu, x('r_1810')/mu, x('r_1672')/glc];
        end
        T = [T; table(v, strains(k,1), {strjoin(ko(inModel), ' ')}, {strjoin(ko(~inModel), ' ')}, ...
            vals(1), vals(2), vals(3), vals(4), vals(5), vals(6), vals(7), vals(8), 'VariableNames', ...
            {'variant','strain','genes_in_model','genes_not_in_model','mu','Y_biomass','Y_ethanol','Y_glycerol', ...
             'acetate_per_biomass','succinate_per_biomass','glycine_per_biomass','CO2_per_glucose'})]; %#ok<AGROW>
    end
end
writetable(T, fullfile(paths.results, outName), 'FileType','text','Delimiter','\t');
disp(T(strcmp(T.variant, variants{1}), :));
end

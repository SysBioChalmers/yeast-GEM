function whatIfC1Evidence()
% whatIfC1Evidence  Bound-free curation candidates for the glycine surplus
%   (no change to the released model). Glycine cleavage (GCV) is reopened in
%   every package; the packages differ only in reaction directions and
%   cofactors, never in flux values:
%     G   GCV open
%     G+H     Hom6 homoserine dehydrogenase NADPH only (r_0546, NADH, removed)
%     G+M     mitochondrial C1-THF dehydrogenase (Mis1, r_0733) reversible, so
%             mitochondrial CH2-THF can be oxidised to formate
%     G+M+S   G+M, mitochondrial SHMT (Shm1, r_0503) serine -> glycine only
%     G+M+S+H all of the above
%     G+T     threonine aldolase (Gly1, r_1040) blocked (reference only)
%     G+H+F   G+H, FMNH2-linked fumarate reductase (r_1000) removed
%     G+H+L   G+H, NAD(P)H:FMN reductase (Lot6, r_0441/r_0442) removed
%     released+F  released model, r_1000 removed
%     released+L+E  released model, Lot6 free-FMN reduction removed, Iah1 esterases hydrolysis only
%   Scored on the wild type, adh3 and gpd1 gpd2, the forced glycine export
%   at optimal growth, the fluxomics benchmark and the growth benchmark.

paths = setupStudy('3.0.0-rc1');
out = fullfile(paths.results, 'whatif_c1');
if ~isfolder(out), mkdir(out); end
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
base = loadTaggedModel(paths.yeast, '9.1.0');
gcv = {'r_0501','r_0507','r_0509'};
G = setParam(base, 'ub', gcv, 1000);
H = @(m) setParam(m, 'ub', 'r_0546', 0);
M = @(m) setParam(m, 'lb', 'r_0733', -1000);
S = @(m) setParam(m, 'lb', 'r_0503', 0);
F = @(m) setParam(setParam(m, 'lb', 'r_1000', 0), 'ub', 'r_1000', 0);   % FMN-linked fumarate reductase removed
L = @(m) setParam(m, 'ub', {'r_0441','r_0442'}, 0);                     % free FMN reduction by Lot6 removed
E = @(m) setParam(m, 'ub', {'r_4711','r_4712','r_4715','r_4716','r_4717','r_4720','r_4721'}, 0);  % no ester synthesis by Iah1
pk = {'v9.1.0 as released', base; 'G', G; 'G+H', H(G); 'G+M', M(G); 'G+M+S', S(M(G)); ...
      'G+M+S+H', H(S(M(G))); 'G+T', setParam(G, 'ub', 'r_1040', 0); ...
      'G+H+F', F(H(G)); 'G+H+L', L(H(G)); 'released+F', F(base); 'released+L+E', E(L(base))};

T = table();
for p = 1:size(pk,1)
    an = anaerobicModel(pk{p,2});
    row = table(pk(p,1), 'VariableNames', {'package'});
    for sc = {'sjoberg', {}, 23; 'bakker_wt', {}, 19.1; 'bakker_adh3', {'YMR083W'}, 12.6; 'gpd1gpd2', {'YDL022W','YOL059W'}, 19.1}'
        mm = an; if ~isempty(sc{2}), mm = removeGenes(mm, sc{2}); end
        mm = setParam(setParam(mm, 'eq', 'r_1714', -sc{3}), 'obj', 'r_4041', 1);
        s = solveLP(mm, 1);
        if isempty(s.x), g = 0; gl = NaN; else
            g = s.x(strcmp(mm.rxns,'r_4041')); gl = s.x(strcmp(mm.rxns,'r_1808')); end
        row.([sc{1} '_growth']) = g; row.([sc{1} '_glycerol']) = gl;
        if strcmp(sc{1}, 'sjoberg')
            row.gcv = s.x(strcmp(mm.rxns,'r_0501')); row.gly1 = s.x(strcmp(mm.rxns,'r_1040'));
            row.shm1 = s.x(strcmp(mm.rxns,'r_0503')); row.mis1_dh = s.x(strcmp(mm.rxns,'r_0733'));
            row.formate_mito_out = s.x(strcmp(mm.rxns,'r_1794'));
            row.succinate = s.x(strcmp(mm.rxns,'r_2056')); row.frd_fmn = s.x(strcmp(mm.rxns,'r_1000'));
            mf = setParam(setParam(mm, 'lb', 'r_4041', 0.99999*g), 'obj', 'r_1810', -1);
            sf = solveLP(mf); row.gly_export_min = sf.x(strcmp(mf.rxns,'r_1810'));
        end
        if strcmp(sc{1}, 'gpd1gpd2') && g > 1e-6
            [~, exch] = getExchangeRxns(mm);
            big = exch(s.x(exch) > 0.05 & ~ismember(mm.rxns(exch), {'r_1761','r_1672','r_2100','r_1832'}));
            row.gpd_exports = {strjoin(strcat(mm.rxnNames(big)', '=', arrayfun(@(v) sprintf('%.2f',v), s.x(big)', 'UniformOutput', false)), ', ')};
        else
            row.gpd_exports = {''};
        end
    end
    figure('Visible','off'); [R2f, mref] = fluxPredSim(an); close all
    figure('Visible','off'); [~, e, gg] = growthNewSim(pk{p,2}); close all
    o = double(e(:,4)); pr = gg(:,4);
    row.flux_R2 = R2f; row.flux_MRE = mref;
    row.growth_R2 = 1 - sum((o-pr).^2)/sum((o-mean(o)).^2); row.growth_MRE = mean(abs(pr-o)./o);
    grp = [repmat({'Nlim_aer'},9,1); repmat({'Clim_aer'},11,1); repmat({'Clim_anaer'},6,1); repmat({'Nlim_anaer'},6,1)];
    for gname = unique(grp)'
        sel = strcmp(grp, gname{1});
        row.(['MRE_' gname{1}]) = mean(abs(pr(sel)-o(sel))./o(sel));
    end
    T = [T; row]; %#ok<AGROW>
end
writetable(T, fullfile(out, 'packages_evidence.tsv'), 'FileType','text','Delimiter','\t');
disp(T(:, {'package','sjoberg_growth','sjoberg_glycerol','gly_export_min','gcv','gly1','shm1','mis1_dh','formate_mito_out','succinate','frd_fmn'}));
disp(T(:, {'package','bakker_wt_growth','bakker_wt_glycerol','bakker_adh3_growth','bakker_adh3_glycerol','gpd1gpd2_growth'}));
disp(T(:, {'package','flux_R2','flux_MRE','growth_R2','growth_MRE','MRE_Clim_aer','MRE_Nlim_aer','MRE_Clim_anaer','MRE_Nlim_anaer'}));
disp(T(:, {'package','gpd_exports'}));
end

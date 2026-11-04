function whatIfC1()
% whatIfC1  What-if curation packages for the glycine surplus (no change to the
%   released model). Glycine cleavage (GCV; r_0501, r_0507, r_0509) is blocked
%   in v9.1.0, so surplus glycine from serine hydroxymethyltransferase can only
%   leave the cell. Packages, applied to the aerobic model (so they also reach
%   the growth benchmark) before the anaerobic constraints:
%     A  GCV reopened
%     B  GCV limited (sum of the three reactions <= gcvCap); glycine export closed
%     C  B + every other non-product export closed, except glycolaldehyde (a
%        dead-end cofactor byproduct), capped at 0.001
%   Each is scored on the wild-type predictions, gpd1 gpd2 viability, the
%   fluxomics benchmark and the chemostat growth benchmark.

paths = setupStudy('3.0.0-rc1');
out = fullfile(paths.results, 'whatif_c1');
if ~isfolder(out), mkdir(out); end
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
base = loadTaggedModel(paths.yeast, '9.1.0');
[~, exch] = getExchangeRxns(base);
gcv = {'r_0501','r_0507','r_0509'};
products = {'r_1761','r_1672','r_1808','r_1634','r_2056','r_2033','r_1546','r_1551', ...
            'r_2100','r_1832','r_4041','r_2111'};
gcvCap = 0.05;

pk = struct('name', {}, 'model', {});
pk(1).name = 'v9.1.0 as released'; pk(1).model = base;
m = base; for r = gcv, m = setParam(m, 'ub', r{1}, 1000); end
pk(2).name = 'A: GCV reopened'; pk(2).model = m;
m = limitGCV(base, gcv, gcvCap); m = setParam(m, 'ub', 'r_1810', 0);
pk(3).name = sprintf('B: GCV <= %g, glycine export closed', gcvCap); pk(3).model = m;
closeIdx = exch(~ismember(base.rxns(exch), [products {'r_1814'}]));
m.ub(closeIdx) = min(m.ub(closeIdx), 0);
m = setParam(m, 'ub', 'r_1814', 0.001);
pk(4).name = 'C: B + other non-product exports closed'; pk(4).model = m;

T = table();
for p = 1:numel(pk)
    an = anaerobicModel(pk(p).model);
    row = table({pk(p).name}, 'VariableNames', {'package'});
    for sc = {'sjoberg', {}, 23; 'bakker_wt', {}, 19.1; 'bakker_adh3', {'YMR083W'}, 12.6; 'gpd1gpd2', {'YDL022W','YOL059W'}, 19.1}'
        mm = an; if ~isempty(sc{2}), mm = removeGenes(mm, sc{2}); end
        mm = setParam(setParam(mm, 'eq', 'r_1714', -sc{3}), 'obj', 'r_4041', 1);
        s = solveLP(mm, 1);
        if isempty(s.x), g = 0; gl = NaN; gcvF = NaN; else
            g = s.x(strcmp(mm.rxns,'r_4041')); gl = s.x(strcmp(mm.rxns,'r_1808'));
            gcvF = sum(s.x(ismember(mm.rxns, gcv))); end
        row.([sc{1} '_growth']) = g; row.([sc{1} '_glycerol']) = gl; row.([sc{1} '_gcv']) = gcvF;
    end
    figure('Visible','off'); [R2f, mref] = fluxPredSim(an); close all
    figure('Visible','off'); [~, e, g] = growthNewSim(pk(p).model); close all
    o = double(e(:,4)); pr = g(:,4);
    row.flux_R2 = R2f; row.flux_MRE = mref;
    row.growth_R2 = 1 - sum((o-pr).^2)/sum((o-mean(o)).^2); row.growth_MRE = mean(abs(pr-o)./o);
    grp = [repmat({'Nlim_aer'},9,1); repmat({'Clim_aer'},11,1); repmat({'Clim_anaer'},6,1); repmat({'Nlim_anaer'},6,1)];
    for gname = unique(grp)'
        sel = strcmp(grp, gname{1});
        row.(['MRE_' gname{1}]) = mean(abs(pr(sel)-o(sel))./o(sel));
    end
    T = [T; row]; %#ok<AGROW>
end
writetable(T, fullfile(out, 'packages.tsv'), 'FileType','text','Delimiter','\t');
disp(T(:, {'package','sjoberg_growth','sjoberg_glycerol','sjoberg_gcv','bakker_wt_growth','bakker_wt_glycerol','bakker_adh3_growth','bakker_adh3_glycerol'}));
disp(T(:, {'package','gpd1gpd2_growth','flux_R2','flux_MRE','growth_R2','growth_MRE','MRE_Clim_aer','MRE_Nlim_aer','MRE_Clim_anaer','MRE_Nlim_anaer'}));
end

function m = limitGCV(m, gcv, cap)
% Open the three GCV reactions and bound their summed flux: a new metabolite
% is produced by each and drained by one capped exchange.
for r = gcv, m = setParam(m, 'ub', r{1}, 1000); end
m.mets{end+1} = 'gcv_budget'; m.metNames{end+1} = 'GCV flux budget';
m.metComps(end+1) = 1; m.S(end+1, :) = 0;
if isfield(m,'b'), m.b(end+1,:) = 0; end
for f = {'metFormulas','metMiriams','metCharges','metSmiles','metDeltaG'}
    if isfield(m, f{1}), if iscell(m.(f{1})), m.(f{1}){end+1} = []; else, m.(f{1})(end+1) = NaN; end, end
end
m.S(end, ismember(m.rxns, gcv)) = 1;
rxn.rxns = {'gcv_budget_out'}; rxn.mets = {'gcv_budget'}; rxn.stoichCoeffs = -1;
rxn.lb = 0; rxn.ub = cap; rxn.rxnNames = {'GCV flux budget'};
m = addRxns(m, rxn, 1, [], false);
end

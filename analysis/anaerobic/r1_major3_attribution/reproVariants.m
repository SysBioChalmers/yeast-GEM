function reproVariants(variants, outName)
% reproVariants  The manuscript's pFBA results (Figure 1, Figure 3A-D) for any
%   list of model variants (see variantModel), e.g.
%   reproVariants({'v911_new','v911_old','v911c_new','v902_new','v902_old'}, 'repro_v911')
%   Sjoberg batch simulation, anaerobic fluxomics benchmark and chemostat
%   growth benchmark with the v9.1.0 test code; writes to results/<outName>.

paths = setupStudy('3.0.0-rc1');
out = fullfile(paths.results, outName);
if ~isfolder(out), mkdir(out); end

% anaerobicModel (original constraints) reads its data relative to code/otherChanges, as in
% anaerobiosis.m.
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
sjoberg = table(); fluxPoints = table(); fluxSummary = table();
growthPoints = table(); growthSummary = table();
growthGroups = [repmat({'N-limited aerobic'},9,1); repmat({'C-limited aerobic'},11,1);
                repmat({'C-limited anaerobic'},6,1); repmat({'N-limited anaerobic'},6,1)];

for k = 1:numel(variants)
    name = variants{k};
    [modelAn, model, ~, post] = variantModel(paths, name);
    if contains(name, '_new'), growthSim = @growthNewSim; else, growthSim = @growthOldSim; end

    % Figure 1: Sjoberg batch, glucose fixed at 23 mmol/gDW/h
    [sim, data, err] = sjobergSim(modelAn);
    sjoberg = [sjoberg; table(repmat({name},4,1), {'glycerol';'ethanol';'CO2';'growth'}, ...
        data(:), err(:), sim(:), 'VariableNames', {'variant','quantity','measured','sd','simulated'})]; %#ok<AGROW>

    % Figure 3C-D: anaerobic fluxomics, flux ratios relative to glucose
    figure('Visible','off');
    [R2, mre, d, s, n, ds] = fluxPredSim(modelAn); close all
    fluxPoints = [fluxPoints; table(repmat({name},numel(d),1), ds(:), n(:), d(:), s(:), ...
        'VariableNames', {'variant','dataset','rxn','measured','simulated'})]; %#ok<AGROW>
    low = d < 30;
    fluxSummary = [fluxSummary; table({name}, R2, corr(d(:),s(:))^2, mre, ...
        r2id(d(low),s(low)), r2id(d(~low),s(~low)), ...
        'VariableNames', {'variant','R2_identity','R2_pearson','MRE','R2_below30','R2_above30'})]; %#ok<AGROW>

    % Figure 3A-B: chemostat growth benchmark (the growth test applies the
    % matching anaerobic function itself)
    figure('Visible','off');
    [~, e, m] = growthSim(model, false, post); close all
    % The old growth test reports a Pearson R2; use the identity-line R2 for
    % every variant, as the manuscript does.
    R2g = r2id(double(e(:,4)), m(:,4));
    growthPoints = [growthPoints; table(repmat({name},size(e,1),1), growthGroups, double(e(:,4)), m(:,4), ...
        'VariableNames', {'variant','condition','measured','simulated'})]; %#ok<AGROW>
    rel = abs(m(:,4)-double(e(:,4)))./double(e(:,4));
    row = table({name}, R2g, mean(rel), 'VariableNames', {'variant','R2_identity','MRE'});
    for g = unique(growthGroups)'
        row.(['MRE_' regexprep(g{1},'[^A-Za-z]','_')]) = mean(rel(strcmp(growthGroups,g{1})));
    end
    growthSummary = [growthSummary; row]; %#ok<AGROW>
    fprintf('%s done\n', name);
end

writetable(sjoberg,       fullfile(out,'fig1_sjoberg.tsv'),      'FileType','text','Delimiter','\t');
writetable(fluxPoints,    fullfile(out,'fig3cd_points.tsv'),     'FileType','text','Delimiter','\t');
writetable(fluxSummary,   fullfile(out,'fig3cd_summary.tsv'),    'FileType','text','Delimiter','\t');
writetable(growthPoints,  fullfile(out,'fig3ab_points.tsv'),     'FileType','text','Delimiter','\t');
writetable(growthSummary, fullfile(out,'fig3ab_summary.tsv'),    'FileType','text','Delimiter','\t');
disp(sjoberg); disp(fluxSummary); disp(growthSummary);
end

function r = r2id(obs, pred)
% Coefficient of determination about the identity line, as in the v9.1.0 tests.
r = 1 - sum((obs-pred).^2)/sum((obs-mean(obs)).^2);
end

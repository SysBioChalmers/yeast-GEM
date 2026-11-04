function makeFigures(optFrac)
% makeFigures  Figures for the revision from the sampling analysis tables.
%   makeFigures(0.99999)           v9.1.0 and v9.0.2, updated and original
%                                  constraints
%   Fig R1  (a) sampled 95% intervals of the reactions discussed in the text,
%               first variant, Sjoberg batch, against pFBA;
%           (b) fluxomics benchmark scored on the samples (CRPS) per model.
%   Fig R2  independent knockout data (Bakker 2000): growth and glycerol,
%           measured against each model.
%   Written to results/figures as PDF and PNG.

if nargin < 1, optFrac = 0.99999; end
paths = setupStudy('3.0.0-rc1');
out = fullfile(paths.results, 'figures');
if ~isfolder(out), mkdir(out); end
slot = [0.165 0.471 0.839; 0.922 0.408 0.204; 0.106 0.686 0.478; 0.929 0.631 0; 0.910 0.482 0.643];  % validated slots 1-5
in = fullfile(paths.results, sprintf('analysis_opt%g', optFrac));
variants = {'v910_new','v910_old','v902_new','v902_old'};
labels   = {'9.1.0, updated constraints','9.1.0, original constraints','9.0.2, updated constraints','9.0.2, original constraints'};
col = slot([1 2 3 4], :);
release = 'v9.1.0';
ink = [0.20 0.20 0.20]; muted = [0.45 0.45 0.45];
rd = @(f) readtable(fullfile(in, f), 'FileType','text','Delimiter','\t');

%% Figure R1
C = rd('claims.tsv');
C = C(strcmp(C.variant, variants{1}), :);
f = figure('Visible','off','Units','centimeters','Position',[0 0 22 11],'Renderer','painters');
t = tiledlayout(f, 1, 3, 'TileSpacing','loose', 'Padding','compact');
ax = nexttile(t, [1 2]); hold(ax,'on');
n = height(C);
for i = 1:n
    y = n - i + 1;
    plot(ax, [C.q2_5(i) C.q97_5(i)], [y y], '-', 'Color', col(1,:), 'LineWidth', 4);
    plot(ax, C.pfba(i), y, 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', ink, 'LineWidth', 1.2);
end
set(ax, 'YTick', 1:n, 'YTickLabel', flip(C.description), 'TickLabelInterpreter','none', 'FontSize', 7, ...
    'XColor', muted, 'YColor', ink, 'Box','off');
xline(ax, 0, ':', 'Color', muted);
lo = min([C.q2_5; C.pfba]); hi = max([C.q97_5; C.pfba]); pad = 0.06*(hi-lo);
xlim(ax, [lo-pad hi+pad]); ylim(ax, [0.4 n+0.6]);
xlabel(ax, 'flux (mmol gDW^{-1} h^{-1})');
title(ax, ['a   ' release ', optimal solutions (Sjoberg batch)'], 'FontWeight','normal', 'FontSize', 9, 'Color', ink);
legend(ax, {'sampled 95% interval','pFBA'}, 'Location','southoutside', 'Box','off', 'FontSize', 7, 'Orientation','horizontal');

ax = nexttile(t); hold(ax,'on');
chainFile = fullfile(paths.results, 'chains', sprintf('crps_opt%g.tsv', optFrac));
crps = zeros(numel(variants), 1);
if isfile(chainFile)
    % four independent chains: data-set means weighted by the number of
    % compared fluxes (glucose excluded); pooled bar, range over chains in the label
    opts = detectImportOptions(chainFile, 'FileType','text', 'Delimiter','\t');
    opts = setvartype(opts, {'variant','dataset','chain'}, 'char');
    X = readtable(chainFile, opts);
    chains = {'1','2','3','4','pooled'};
    spread = zeros(numel(variants), 1);
    for k = 1:numel(variants)
        per = zeros(1, numel(chains));
        for c = 1:numel(chains)
            r = strcmp(X.variant, variants{k}) & strcmp(X.chain, chains{c});
            w = X.n(r) - 1;
            per(c) = sum(w .* X.crps_mean_no_glucose(r)) / sum(w);
        end
        crps(k) = per(end); spread(k) = max(per(1:end-1)) - min(per(1:end-1));
    end
    xl = {'mean CRPS', '(% of glucose uptake)', sprintf('4 chains, range \\leq %.2f', ceil(100*max(spread))/100)};
else
    S = rd('fluxomics_summary.tsv');
    S = S(~S.with_glucose_points, :);
    for k = 1:numel(variants)
        crps(k) = S.crps_mean(strcmp(S.variant, variants{k}));
    end
    xl = {'mean CRPS','(% of glucose uptake)'};
end
for k = 1:numel(variants)
    y = numel(variants) - k + 1;
    barh(ax, y, crps(k), 0.6, 'FaceColor', col(k,:), 'EdgeColor', 'none');
    text(ax, crps(k), y, sprintf('  %.2f', crps(k)), 'HorizontalAlignment','left', ...
        'VerticalAlignment','middle', 'FontSize', 7, 'Color', ink);
end
set(ax, 'YTick', 1:numel(variants), 'YTickLabel', flip(labels), 'FontSize', 7, 'XColor', muted, 'YColor', ink, 'Box','off');
xlim(ax, [0 1.2*max(crps)]);
xlabel(ax, xl);
title(ax, 'b   Fluxomics score', 'FontWeight','normal', 'FontSize', 9, 'Color', ink);
saveFig(f, fullfile(out, 'figR1_sampling'));

%% Figure R2
K = rd('knockouts.tsv');
f = figure('Visible','off','Units','centimeters','Position',[0 0 14 7.5],'Renderer','painters');
t = tiledlayout(f, 1, 2, 'TileSpacing','loose', 'Padding','compact');
panels = {'mu','growth rate (h^{-1})','a   Growth'; 'q_glycerol','glycerol (mmol gDW^{-1} h^{-1})','b   Glycerol'};
strains = {'bakker_wt','wild type'; 'bakker_adh3','adh3\Delta'};
hv = gobjects(1, numel(variants)); hm = gobjects(1);
for p = 1:2
    ax = nexttile(t); hold(ax,'on');
    for s = 1:2
        rows = strcmp(K.scenario, strains{s,1}) & strcmp(K.quantity, panels{p,1});
        R = K(rows,:);
        m = R.measured(1); sd = R.measured_sd(1);
        patch(ax, [s-0.4 s+0.4 s+0.4 s-0.4], [m-sd m-sd m+sd m+sd], [0.85 0.85 0.85], 'EdgeColor','none');
        hm = plot(ax, [s-0.4 s+0.4], [m m], '-', 'Color', ink, 'LineWidth', 1.5);
        for k = 1:numel(variants)
            v = R.median(strcmp(R.variant, variants{k}));
            hv(k) = plot(ax, s + 0.64*((k-1)/(numel(variants)-1) - 0.5), v, 'o', 'MarkerSize', 7, 'MarkerFaceColor', col(k,:), 'MarkerEdgeColor', 'w', 'LineWidth', 1);
        end
    end
    set(ax, 'XTick', 1:2, 'XTickLabel', strains(:,2), 'XLim', [0.4 2.6], 'FontSize', 8, 'XColor', muted, 'YColor', muted, 'Box','off');
    yl = ylim(ax); ylim(ax, [yl(1) - 0.06*diff(yl), yl(2) + 0.04*diff(yl)]);
    ylabel(ax, panels{p,2});
    title(ax, panels{p,3}, 'FontWeight','normal', 'FontSize', 9, 'Color', ink);
end
lg = legend(hv(1).Parent, [hv hm], [labels, {'measured \pm SD'}], 'Box','off','FontSize',7, 'Orientation','horizontal', 'NumColumns', 2);
lg.Layout.Tile = 'south';
saveFig(f, fullfile(out, 'figR2_knockouts'));
end

function saveFig(f, base)
% Headless-safe export: vector PDF and 300 dpi PNG with the painters renderer.
set(f, 'PaperUnits', 'centimeters');
pos = get(f, 'Position');
set(f, 'PaperSize', pos(3:4), 'PaperPosition', [0 0 pos(3:4)]);
print(f, [base '.pdf'], '-dpdf', '-painters');
print(f, [base '.png'], '-dpng', '-r300', '-painters');
end

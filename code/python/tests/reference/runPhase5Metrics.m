function runPhase5Metrics(yeastGemPath, outJson)
% runPhase5Metrics  Compute the Tier-3 metrics with MATLAB for cross-
% language verification.
%
%   Writes a JSON file with growth R², essential-gene accuracy /
%   sensitivity / specificity / MCC, and anaerobic-flux R² + MRE.
%   The Python equivalent metrics live under yeastgem.model_tests
%   and should match within float tolerance.

warning('off','all');
restoredefaultpath; rehash toolboxcache;
addpath('/opt/gurobi1301/linux64/matlab');
addpath(genpath('/home/eduardk/github/RAVEN'));
addpath(fullfile(yeastGemPath, 'code'));
addpath(fullfile(yeastGemPath, 'code', 'modelCuration'));
addpath(fullfile(yeastGemPath, 'code', 'otherChanges'));
addpath(fullfile(yeastGemPath, 'code', 'missingFields'));
addpath(fullfile(yeastGemPath, 'code', 'modelTests'));

model = loadYeastModel;

% --- growth -----------------------------------------------------------
fprintf('Running growth...\n');
fig = figure('Visible','off');
growth_r2 = growth(model);
close(fig);

% --- essential genes ---------------------------------------------------
fprintf('Running essentialGenes...\n');
[acc, tp, tn, fn, fp] = essentialGenes(model);
n_tp = numel(tp); n_tn = numel(tn); n_fp = numel(fp); n_fn = numel(fn);
sens = 100*n_tp / (n_tp + n_fn);
spec = 100*n_tn / (n_tn + n_fp);
denom_mcc = (n_tp + n_fp)*(n_tp + n_fn)*(n_tn + n_fp)*(n_tn + n_fn);
mcc = (n_tp*n_tn - n_fp*n_fn) / sqrt(denom_mcc);

% --- anaerobic flux ----------------------------------------------------
fprintf('Running anaerobic_flux_predictions...\n');
model_an = applyYeastCondition(model, 'anaerobic');
fig = figure('Visible','off');
prev = cd(fullfile(yeastGemPath, 'code', 'modelTests'));
cleanup = onCleanup(@() cd(prev));
anaerobic_r2 = anaerobic_flux_predictions(model_an);
clear cleanup;
close(fig);

% --- emit JSON --------------------------------------------------------
out = struct();
out.growth_r2 = growth_r2;
out.essential_genes_accuracy = acc;
out.essential_genes_sensitivity = sens;
out.essential_genes_specificity = spec;
out.essential_genes_mcc = mcc;
out.essential_genes_tp = n_tp;
out.essential_genes_tn = n_tn;
out.essential_genes_fp = n_fp;
out.essential_genes_fn = n_fn;
out.anaerobic_flux_r2 = anaerobic_r2;

fid = fopen(outJson, 'w');
fprintf(fid, '%s', jsonencode(out, 'PrettyPrint', true));
fclose(fid);
fprintf('Wrote %s\n', outJson);
end

function curationAblation(outName, fine)
% curationAblation  Which part of the 9.0.2 -> 9.1.0 curation changes the
%   anaerobic predictions? Replays code/modelCuration/v9_1_0.m section by
%   section on yeast-GEM 9.0.2: cumulatively (sections 1..k) and leaving one
%   section out. Every model gets the updated anaerobic constraints (9.1.0
%   anaerobicModel) and the manuscript benchmarks: Sjoberg batch (Figure 1),
%   anaerobic fluxomics (Figure 3C-D) and chemostat growth (Figure 3A-B).
%   Writes results/<outName>.tsv.

if nargin < 1, outName = 'curation_ablation'; end
if nargin < 2, fine = false; end
paths = setupStudy('3.0.0-rc1');
scriptFile = fullfile(paths.yeast,'code','modelCuration','v9_1_0.m');
L = splitlines(fileread(scriptFile));
% Sections of v9_1_0.m (line ranges, inclusive)
S = {
    'Oac1 malate/sulphate transport',           [16 51]
    'H+ symporters irreversible',               [52 70]
    'mass and charge balance',                  [71 174]
    'Gcy1 reversed and irreversible',           [175 186]
    'mitochondrial ATP synthase reversible',    [187 191]
    'C1/folate directions',                     [195 205]
    'other directions and cofactors',           [192 194; 206 222]
    'balancing and removals',                   [223 319]
    'GCV, Glt1, carnitine blocked',             [320 336]
    'amino-acid ratio (Bjorkeroth 2020)',       [337 340]
    'protein rescaled to 1 g/gDW',              [341 346]
    'biomass pseudoreaction H+ balance',        [347 370]
    'biomass NADPH (degree of reduction)',      [371 381]
    'glycine export allowed',                   [382 388]
    };
if fine   % split sections 7 and 9; leave-one-out only
    S = [S(1:6,:);
        {'r_0227 renamed', [192 194]
         'PFK1 side reaction r_0887 removed', [206 208]
         'TYR1 cofactor (r_0939)', [209 214]
         'esterification irreversible', [215 218]
         'polyphosphate irreversible', [219 222]};
        S(8,:);
        {'GCV blocked', [320 327]
         'Glt1 blocked', [328 331]
         'carnitine acetyltransferase blocked', [332 336]};
        S(10:end,:)];
end
code = cell(size(S,1),1);
for i = 1:size(S,1)
    r = S{i,2}; c = {};
    for j = 1:size(r,1), c = [c; L(r(j,1):r(j,2))]; end %#ok<AGROW>
    code{i} = strjoin(c, newline);
end
assert(contains(L{389},'%% ====') && contains(L{16},'%% ===='), 'v9_1_0.m line layout changed');

base902 = loadTaggedModel(paths.yeast,'9.0.2');
ref910  = loadTaggedModel(paths.yeast,'9.1.0');
n = size(S,1);
sets = {'9.0.2 (no section)', []};
if ~fine
    for k = 1:n, sets(end+1,:) = {sprintf('cumulative 1-%d (+ %s)',k,S{k,1}), 1:k}; end %#ok<AGROW>
end
sets(end+1,:) = {'all sections', 1:n};
for k = 1:n, sets(end+1,:) = {sprintf('all but %d (%s)',k,S{k,1}), setdiff(1:n,k)}; end %#ok<AGROW>
sets(end+1,:) = {'released 9.1.0', 'ref'};

T = table();
for s = 1:size(sets,1)
    name = sets{s,1};
    if ischar(sets{s,2})
        model = ref910; ok = true; msg = '';
    else
        [model, ok, msg] = replay(base902, code(sets{s,2}), paths);
    end
    row = table({name}, ok, {msg}, 'VariableNames', {'set','ok','error'});
    vals = nan(1,10);
    if ok
        try
            vals = benchmarks(model, paths);
        catch e
            row.ok = false; row.error = {e.message};
        end
    end
    row = [row, array2table(vals, 'VariableNames', {'growth','glycerol','ethanol','CO2', ...
        'flux_R2','flux_MRE','growth_R2','growth_MRE','MRE_Clim_anaerobic','MRE_Nlim_anaerobic'})]; %#ok<AGROW>
    T = [T; row]; %#ok<AGROW>
    fprintf('%-60s glycerol %.3f  flux MRE %.4f  growth R2 %.4f\n', name, vals(2), vals(6), vals(7));
end
writetable(T, fullfile(paths.results,[outName '.tsv']), 'FileType','text','Delimiter','\t');
end

function [model, ok, msg] = replay(model, sections, paths)
% Runs v9_1_0.m sections in code/modelCuration, as the script does.
ok = true; msg = '';
here = cd(fullfile(paths.yeast,'code','modelCuration'));
back = onCleanup(@() cd(here));
codeDir = fullfile(paths.yeast,'code'); %#ok<NASGU>
for i = 1:numel(sections)
    try
        evalc(sections{i});
    catch e
        ok = false; msg = e.message; return
    end
end
end

function vals = benchmarks(model, paths)
here = cd(fullfile(paths.yeast,'code','otherChanges'));
back = onCleanup(@() cd(here));
modelAn = anaerobicModel(model);
sim = sjobergSim(modelAn);
figure('Visible','off'); [R2f, mref] = fluxPredSim(modelAn); close all
figure('Visible','off'); [~, e, m] = growthNewSim(model, false); close all
obs = double(e(:,4)); pred = m(:,4);
R2g = 1 - sum((obs-pred).^2)/sum((obs-mean(obs)).^2);
rel = abs(pred-obs)./obs;
vals = [sim(4) sim(1) sim(2) sim(3) R2f mref R2g mean(rel) mean(rel(21:26)) mean(rel(27:32))];
end

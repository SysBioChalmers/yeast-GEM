function curationAblation(outName, mode, only)
% curationAblation  Which part of the 9.0.2 -> 9.1.0 curation changes the
%   anaerobic predictions? Replays code/modelCuration/v9_1_0.m section by
%   section on yeast-GEM 9.0.2: cumulatively (sections 1..k) and leaving one
%   section out. Every model gets the updated anaerobic constraints (9.1.0
%   applyEnvironment(model,'anaerobic')) and the manuscript benchmarks: Sjoberg batch (Figure 1),
%   anaerobic fluxomics (Figure 3C-D) and chemostat growth (Figure 3A-B).
%   mode 'sections' (default) also gives the cumulative models; 'fine' splits
%   sections 7 and 9 and leaves out one at a time; 'categories' adds the leave-out
%   and add-in of groups of sections (the categories of the reviewer comment) and
%   the changes of 9.1.2 after 9.1.0 (thiamine, fumarate reductases, TYR1 and the acetate transport and ribonucleotide reductase reactions of #389), which
%   are applied from the checkout model on top of the replay.
%   mode 'pairs' leaves out every two sections (190 models) and every two categories (10 models).
%   only: indices of the models to compute (for parallel jobs); the file is then
%   results/<outName>_<index>.tsv. Writes results/<outName>.tsv.

if nargin < 1, outName = 'curation_ablation'; end
if nargin < 2, mode = 'sections'; end
if nargin < 3, only = []; end
fine = ismember(mode, {'fine','categories','pairs'});
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
final = [];
GROUPS = struct('thiamine', {{'r_1603','r_4788'}}, 'fumarate', {{'r_0454','r_0455','r_1000','r_4264'}}, 'tyr1', {{'r_0939'}}, ...
    'acetate_rnr', {{'r_4785','r_4786','r_4787'}});   % last: mitochondrial acetate transport and ribonucleotide reductase (#389)
allGroups = fieldnames(GROUPS)';
cats = {'charge and element balance', [3 12]
        'biomass equation', [16 17 18 19]
        'mitochondrial and nitrogen metabolism', [1 5 13 14 15 20]
        'reaction directions and cofactors', [2 4 6 7 8 9 10 11]};
if ismember(mode, {'categories','pairs'}), final = loadTaggedModel(paths.yeast, 'dev'); end
sets = {'9.0.2 (no section)', [], {}};
if strcmp(mode, 'pairs'), sets = cell(0,3); end
if ~fine
    for k = 1:n, sets(end+1,:) = {sprintf('cumulative 1-%d (+ %s)',k,S{k,1}), 1:k, {}}; end %#ok<AGROW>
end
grp = {}; if ismember(mode, {'categories','pairs'}), grp = allGroups; end
sets(end+1,:) = {'all sections', 1:n, grp};
for k = 1:n, sets(end+1,:) = {sprintf('all but %d (%s)',k,S{k,1}), setdiff(1:n,k), grp}; end %#ok<AGROW>
if strcmp(mode, 'categories')
    assert(n == 20, 'category indices assume the 20 sections of the fine split');
    for c = 1:size(cats,1)
        sets(end+1,:) = {['all but category: ' cats{c,1}], setdiff(1:n, cats{c,2}), grp}; %#ok<AGROW>
        sets(end+1,:) = {['only category: ' cats{c,1}], cats{c,2}, {}}; %#ok<AGROW>
    end
    sets(end+1,:) = {'all but category: changes of 9.1.2 after 9.1.0', 1:n, {}};
    sets(end+1,:) = {'only category: changes of 9.1.2 after 9.1.0', [], allGroups};
    for g = allGroups, sets(end+1,:) = {['all but 9.1.2 change: ' g{1}], 1:n, setdiff(allGroups, g)}; end %#ok<AGROW>
    sets(end+1,:) = {'checkout (9.1.2)', 'dev', {}};
end
if strcmp(mode, 'pairs')   % every two sections, and every two categories, left out together
    assert(n == 20, 'category indices assume the 20 sections of the fine split');
    sets = cell(0,3);
    for i = 1:n-1, for j = i+1:n
        sets(end+1,:) = {sprintf('all but %d and %d (%s; %s)', i, j, S{i,1}, S{j,1}), setdiff(1:n,[i j]), allGroups}; end, end %#ok<AGROW>
    cats(end+1,:) = {'changes of 9.1.2 after 9.1.0', []};
    for a = 1:size(cats,1)-1, for b = a+1:size(cats,1)
        gr = allGroups; if b == size(cats,1), gr = {}; end
        sets(end+1,:) = {['all but categories: ' cats{a,1} '; ' cats{b,1}], setdiff(1:n,[cats{a,2} cats{b,2}]), gr}; end, end %#ok<AGROW>
else
    sets(end+1,:) = {'released 9.1.0', 'ref', {}};
end
if isempty(only), only = 1:size(sets,1); end

T = table();
for s = only
    name = sets{s,1};
    if ischar(sets{s,2})
        if strcmp(sets{s,2}, 'dev'), model = final; else, model = ref910; end
        ok = true; msg = '';
    else
        [model, ok, msg] = replay(base902, code(sets{s,2}), paths);
        for g = sets{s,3}, model = copyReactions(model, final, GROUPS.(g{1})); end
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
if numel(only) < size(sets,1), outName = sprintf('%s_%d', outName, only(1)); end
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
modelAn = applyEnvironment(model, 'anaerobic');
sim = sjobergSim(modelAn);
figure('Visible','off'); [R2f, mref] = fluxPredSim(modelAn); close all
figure('Visible','off'); [~, e, m] = growthNewSim(model, false); close all
obs = double(e(:,4)); pred = m(:,4);
R2g = 1 - sum((obs-pred).^2)/sum((obs-mean(obs)).^2);
rel = abs(pred-obs)./obs;
vals = [sim(4) sim(1) sim(2) sim(3) R2f mref R2g mean(rel) mean(rel(21:26)) mean(rel(27:32))];
end

function model = copyReactions(model, final, ids)
% Reactions of the checkout model (9.1.2) put into the replayed model: the
% stoichiometry and bounds of existing reactions are replaced, new reactions
% and metabolites are added, and reactions that the checkout model lacks are removed.
for id = ids
    fj = find(strcmp(final.rxns, id{1}));
    if isempty(fj)   % removed in the checkout model
        model = removeReactions(model, id{1}, true);
        continue
    end
    fi = find(final.S(:, fj));
    metIds = final.mets(fi); coef = full(final.S(fi, fj));
    missing = ~ismember(metIds, model.mets);
    if any(missing)
        add.mets = metIds(missing); add.metNames = final.metNames(fi(missing));
        add.compartments = final.comps(final.metComps(fi(missing)));
        model = addMets(model, add);
    end
    mj = find(strcmp(model.rxns, id{1}));
    if isempty(mj)
        r.rxns = id; r.rxnNames = final.rxnNames(fj); r.mets = {metIds'}; r.stoichCoeffs = {coef'};
        r.lb = final.lb(fj); r.ub = final.ub(fj);
        model = addRxns(model, r);
    else
        [~, loc] = ismember(metIds, model.mets);
        model.S(:, mj) = 0; model.S(loc, mj) = coef;
        model.lb(mj) = final.lb(fj); model.ub(mj) = final.ub(fj);
    end
end
end

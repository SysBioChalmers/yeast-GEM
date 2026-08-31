function migrateAnnotationTsvs()
% migrateAnnotationTsvs
%   One-time migration for yeast-GEM#379: splits reaction, metabolite and
%   gene cross-reference annotation out of model/yeast-GEM.yml into
%   model/reactions.tsv, model/metabolites.tsv and model/genes.tsv, and
%   rewrites model/yeast-GEM.yml without those fields. sbo, deltaG,
%   confidence_score, notes and subsystem are untouched.
%
%   Verifies itself: reloads the rewritten yml, merges the new tsvs back
%   in via annotateGEM.m, and confirms the result is identical to the
%   original load for every reaction/metabolite/gene, not sampled.
%
%   Kept as a historical record, like the code/modelCuration/vX_Y_Z.m
%   curation scripts -- not meant to be run again.
%
% Usage:
%   migrateAnnotationTsvs()

modelPath = fullfile(fileparts(mfilename('fullpath')),'..','..','model');
ymlFile = fullfile(modelPath,'yeast-GEM.yml');

fprintf('Loading %s\n', ymlFile);
original = readYAMLmodel(ymlFile);

%% Reactions: name (read-only convenience copy, yml stays authoritative --
% see model/README.md) + bigg.reaction / kegg.pathway / kegg.reaction /
% metanetx.reaction (miriams) + ec-code
rxnMiriamCols = {'bigg.reaction','kegg.pathway','kegg.reaction','metanetx.reaction'};
[rxnTab, leanRxnMiriams] = extractAndStripMiriams(original.rxns, original.rxnMiriams, rxnMiriamCols);
rxnTab.name = original.rxnNames;
if isfield(original,'eccodes')
    rxnTab.("ec-code") = original.eccodes;
else
    rxnTab.("ec-code") = repmat({''}, numel(original.rxns), 1);
end
leanEccodes = repmat({''}, numel(original.rxns), 1);

writeAnnotationTsv(fullfile(modelPath,'reactions.tsv'), rxnTab, ...
    {'id','name','bigg.reaction','ec-code','kegg.pathway','kegg.reaction','metanetx.reaction'});

%% Metabolites: name (same read-only-copy convention) + bigg.metabolite /
% chebi / kegg.compound / metanetx.chemical (miriams) + smiles
metMiriamCols = {'bigg.metabolite','chebi','kegg.compound','metanetx.chemical'};
[metTab, leanMetMiriams] = extractAndStripMiriams(original.mets, original.metMiriams, metMiriamCols);
metTab.name = original.metNames;
if isfield(original,'metSmiles')
    metTab.smiles = original.metSmiles;
else
    metTab.smiles = repmat({''}, numel(original.mets), 1);
end
leanSmiles = repmat({''}, numel(original.mets), 1);

writeAnnotationTsv(fullfile(modelPath,'metabolites.tsv'), metTab, ...
    {'id','name','bigg.metabolite','chebi','kegg.compound','metanetx.chemical','smiles'});

%% Genes: uniprot (miriam only)
[geneTab, leanGeneMiriams] = extractAndStripMiriams(original.genes, original.geneMiriams, {'uniprot'});

writeAnnotationTsv(fullfile(modelPath,'genes.tsv'), geneTab, {'id','uniprot'});

%% Verify the name columns just written, read back from disk rather than
% trusting the in-memory assignment above -- catches a writetable/quoting
% surprise the direct assignment wouldn't.
verifyTsvNamesMatch(fullfile(modelPath,'reactions.tsv'), original.rxns, original.rxnNames, 'reaction');
verifyTsvNamesMatch(fullfile(modelPath,'metabolites.tsv'), original.mets, original.metNames, 'metabolite');
fprintf('Name columns verified against model/yeast-GEM.yml.\n');

%% Build & write the lean model
lean = original;
lean.rxnMiriams  = leanRxnMiriams;
lean.eccodes     = leanEccodes;
lean.metMiriams  = leanMetMiriams;
lean.metSmiles   = leanSmiles;
lean.geneMiriams = leanGeneMiriams;

fprintf('Writing lean %s\n', ymlFile);
writeYAMLmodel(lean, ymlFile);

%% Verify: reload the lean yml, merge the tsvs back via annotateGEM, and
% confirm the result is identical to the original load.
fprintf('Verifying round trip...\n');
reloaded = readYAMLmodel(ymlFile);
merged = annotateGEM(reloaded, modelPath);

verifyRoundTrip(original, merged);

fprintf('Migration verified: lossless.\n');
end

function [tab, leanMiriams] = extractAndStripMiriams(ids, miriamCellArray, columns)
% Split miriamCellArray (model.{rxn,met,gene}Miriams) into: a table of
% one column per requested namespace (';'-joined for multi-value
% entries), and a copy of miriamCellArray with those namespaces removed
% (any other namespace, e.g. sbo, is left untouched; an entity left with
% no annotation at all collapses to '', matching readYAMLmodel's own
% empty-entry convention rather than an empty-but-present struct).
n = numel(ids);
tab = table(ids, 'VariableNames', {'id'});
for c = 1:numel(columns)
    tab.(columns{c}) = repmat({''}, n, 1);
end

leanMiriams = miriamCellArray;
for i = 1:n
    if ~isstruct(miriamCellArray{i})
        continue
    end
    names  = miriamCellArray{i}.name;
    values = miriamCellArray{i}.value;
    keep = true(numel(names),1);
    for c = 1:numel(columns)
        isCol = strcmp(names, columns{c});
        if any(isCol)
            tab.(columns{c}){i} = strjoin(values(isCol), ';');
        end
        keep = keep & ~isCol;
    end
    if all(keep)
        leanMiriams{i} = miriamCellArray{i};
    elseif any(keep)
        leanMiriams{i} = struct('name',{names(keep)},'value',{values(keep)});
    else
        leanMiriams{i} = '';
    end
end
end

function writeAnnotationTsv(filename, tab, columnOrder)
tab = tab(:, columnOrder);
writetable(tab, filename, 'FileType','text', 'Delimiter','\t');
fprintf('Wrote %s (%d rows)\n', filename, height(tab));
end

function verifyTsvNamesMatch(tsvFile, ids, names, label)
opts = detectImportOptions(tsvFile, 'FileType','text', 'Delimiter','\t', ...
    'VariableNamingRule','preserve');
opts = setvartype(opts, opts.VariableNames, 'char');
tab = readtable(tsvFile, opts);
[hasEntry, tsvIdx] = ismember(ids, tab.id);
if ~all(hasEntry)
    error('migrateAnnotationTsvs:verify', '%s: %d id(s) missing from %s', ...
        label, sum(~hasEntry), tsvFile);
end
tsvNames = tab.name(tsvIdx);
mismatch = find(~strcmp(names, tsvNames), 1);
if ~isempty(mismatch)
    error('migrateAnnotationTsvs:verify', ...
        '%s %s: name mismatch between model (%s) and %s (%s)', ...
        label, ids{mismatch}, names{mismatch}, tsvFile, tsvNames{mismatch});
end
end

function verifyRoundTrip(original, merged)
checkMiriamsMatch(original.rxns,  original.rxnMiriams,  merged.rxnMiriams,  'reaction');
checkMiriamsMatch(original.mets,  original.metMiriams,  merged.metMiriams,  'metabolite');
checkMiriamsMatch(original.genes, original.geneMiriams, merged.geneMiriams, 'gene');

% isequaln, not isequal: rxnDeltaG/metDeltaG/metCharges/rxnConfidenceScores
% legitimately contain NaN for unknown values, and isequal treats NaN as
% unequal to itself -- confirmed directly (isequal(model.rxnDeltaG,
% model.rxnDeltaG) on an unmodified reload is already false, purely from
% its 121 NaN entries), so isequal here would flag entities that never
% changed at all.
if ~isequaln(original.eccodes, merged.eccodes)
    error('migrateAnnotationTsvs:verify', 'eccodes mismatch after round trip')
end
if ~isequaln(original.metSmiles, merged.metSmiles)
    error('migrateAnnotationTsvs:verify', 'metSmiles mismatch after round trip')
end

% Sanity: nothing outside the six migrated fields moved.
otherFields = {'mets','rxns','genes','rxnNames','metNames','S','lb','ub', ...
    'grRules','subSystems','metFormulas','metCharges','rxnDeltaG','metDeltaG', ...
    'rxnConfidenceScores','rxnNotes','metNotes','comps','compNames'};
for f = otherFields
    if isfield(original, f{1}) && isfield(merged, f{1})
        if ~isequaln(original.(f{1}), merged.(f{1}))
            error('migrateAnnotationTsvs:verify', 'Unexpected difference in field: %s', f{1})
        end
    elseif isfield(original, f{1}) ~= isfield(merged, f{1})
        error('migrateAnnotationTsvs:verify', 'Field presence differs: %s', f{1})
    end
end
end

function checkMiriamsMatch(ids, miriamsA, miriamsB, label)
for i = 1:numel(ids)
    triA = sort(flattenMiriam(miriamsA{i}));
    triB = sort(flattenMiriam(miriamsB{i}));
    if ~isequal(triA, triB)
        error('migrateAnnotationTsvs:verify', '%s %s: miriam mismatch\n  original: %s\n  merged:   %s', ...
            label, ids{i}, strjoin(triA,'; '), strjoin(triB,'; '))
    end
end
end

function tri = flattenMiriam(m)
if ~isstruct(m)
    tri = {};
    return
end
tri = strcat(m.name, {'='}, m.value);
end

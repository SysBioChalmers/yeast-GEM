function model = annotateGEM(model, annPath)
% annotateGEM
%   Merge reaction, metabolite and gene cross-reference annotation from
%   model/reactions.tsv, model/metabolites.tsv and model/genes.tsv into a
%   model loaded from the (annotation-light) yeast-GEM.yml.
%
%   Called by loadYeastYaml.m so that day-to-day curation sees the same
%   fully-annotated model as before the cross-reference annotation was
%   split out of the yml (yeast-GEM#379). Does not touch SBO terms
%   (addSBOterms.m, unrelated to this migration) or any field other than
%   the ones listed below.
%
% Input:
%   model     model structure, as returned by readYAMLmodel.
%   annPath   (opt) folder containing reactions.tsv, metabolites.tsv and
%             genes.tsv. Default: the model/ folder next to this
%             function's repo root.
%
% Output:
%   model     model structure with the following merged in, by id:
%               reactions:   bigg.reaction, kegg.pathway, kegg.reaction,
%                            metanetx.reaction (into rxnMiriams),
%                            ec-code (into eccodes)
%               metabolites: bigg.metabolite, chebi, kegg.compound,
%                            metanetx.chemical (into metMiriams),
%                            smiles (into metSmiles)
%               genes:       uniprot (into geneMiriams)
%
% Usage:
%   model = annotateGEM(model, annPath);

if nargin < 2 || isempty(annPath)
    funcDir = dbstack('-completenames');
    funcDir = regexprep(funcDir(1).file,[funcDir(1).name '\.m'],'');
    annPath = fullfile(funcDir,'..','model');
end

model = mergeMiriams(model, 'rxnMiriams', 'rxns', fullfile(annPath,'reactions.tsv'), ...
    {'bigg.reaction','kegg.pathway','kegg.reaction','metanetx.reaction'});
model = mergeScalarField(model, 'eccodes', 'rxns', fullfile(annPath,'reactions.tsv'), 'ec-code');

model = mergeMiriams(model, 'metMiriams', 'mets', fullfile(annPath,'metabolites.tsv'), ...
    {'bigg.metabolite','chebi','kegg.compound','metanetx.chemical'});
model = mergeScalarField(model, 'metSmiles', 'mets', fullfile(annPath,'metabolites.tsv'), 'smiles');

model = mergeMiriams(model, 'geneMiriams', 'genes', fullfile(annPath,'genes.tsv'), {'uniprot'});

end

function tab = readAnnotationTsv(filename)
% Read a tsv with every column forced to text, so a column that happens
% to look all-numeric (or is entirely empty) is never silently
% auto-typed as double — every value here is an identifier, never a
% number to compute with. VariableNamingRule 'preserve' keeps dotted
% namespace headers (e.g. "kegg.reaction") intact; they're then read via
% dynamic field access (tab.("kegg.reaction")), which — unlike the
% static tab.kegg.reaction syntax — accepts any header text verbatim.
opts = detectImportOptions(filename, 'FileType','text', 'Delimiter','\t', ...
    'VariableNamingRule','preserve');
opts = setvartype(opts, opts.VariableNames, 'char');
tab = readtable(filename, opts);
end

function model = mergeScalarField(model, fieldName, idField, tsvFile, column)
% Merge a single tsv column into a plain per-entity string field
% (eccodes, metSmiles) that already stores its value in exactly the tsv's
% own format (';'-joined, unquoted) -- a direct assignment, no
% (name,value) pair structure involved.
tab = readAnnotationTsv(tsvFile);
[hasEntry, tsvIdx] = ismember(model.(idField), tab.id);
if ~isfield(model, fieldName)
    model.(fieldName) = repmat({''}, numel(model.(idField)), 1);
elseif numel(model.(fieldName)) < numel(model.(idField))
    model.(fieldName)(end+1:numel(model.(idField)),1) = {''};
end
for i = 1:numel(model.(idField))
    if ~hasEntry(i)
        continue
    end
    value = strtrim(tab.(column){tsvIdx(i)});
    if ~isempty(value)
        model.(fieldName){i,1} = value;
    end
end
end

function model = mergeMiriams(model, miriamField, idField, tsvFile, columns)
% Merge one or more tsv columns into a {name,value} MIRIAM struct array,
% appending to whatever is already there (e.g. an existing 'sbo' entry)
% rather than replacing it. Every column value is split on ';', so a
% multi-valued column (e.g. kegg.pathway) becomes one MIRIAM entry per
% value -- the same shape readYAMLmodel.m itself produces, which is what
% writeYAMLmodel.m expects.
tab = readAnnotationTsv(tsvFile);
[hasEntry, tsvIdx] = ismember(model.(idField), tab.id);
if ~isfield(model, miriamField)
    model.(miriamField) = repmat({''}, numel(model.(idField)), 1);
elseif numel(model.(miriamField)) < numel(model.(idField))
    model.(miriamField)(end+1:numel(model.(idField)),1) = {''};
end

for i = 1:numel(model.(idField))
    if ~hasEntry(i)
        continue
    end
    row = tsvIdx(i);
    if isstruct(model.(miriamField){i})
        names  = model.(miriamField){i}.name;
        values = model.(miriamField){i}.value;
    else
        names  = cell(0,1);
        values = cell(0,1);
    end
    for c = 1:numel(columns)
        cellValue = strtrim(tab.(columns{c}){row});
        if isempty(cellValue)
            continue
        end
        parts = strtrim(strsplit(cellValue, ';'));
        parts = parts(~cellfun(@isempty, parts));
        for p = 1:numel(parts)
            names{end+1,1}  = columns{c};
            values{end+1,1} = parts{p};
        end
    end
    if ~isempty(names)
        model.(miriamField){i,1} = struct('name',{names},'value',{values});
    end
end
end

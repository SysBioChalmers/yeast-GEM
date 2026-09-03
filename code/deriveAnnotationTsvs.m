function model = deriveAnnotationTsvs(model, annPath)
% deriveAnnotationTsvs
%   Write model/reactions.tsv, model/metabolites.tsv and model/genes.tsv
%   from the model's current annotation (rxnMiriams/eccodes, metMiriams/
%   metSmiles, geneMiriams) -- the reverse of annotateGEM. Each tsv is
%   rewritten from scratch, one row per id, in the model's current order.
%
%   Closes a gap where a cross-reference edit made programmatically (e.g.
%   via editMiriam), rather than by hand-editing a tsv cell, would
%   otherwise be silently discarded the next time the model is saved
%   (stripTsvAnnotation strips these fields from model/yeast-GEM.yml and
%   never writes them anywhere else).
%
% Input:
%   model     model structure, with rxnMiriams/metMiriams/geneMiriams
%             and related fields as merged in by annotateGEM.
%   annPath   (opt) folder to write reactions.tsv, metabolites.tsv and
%             genes.tsv into. Default: the model/ folder next to this
%             function's repo root.
%
% Usage:
%   model = deriveAnnotationTsvs(model, annPath);

if nargin < 2 || isempty(annPath)
    funcDir = dbstack('-completenames');
    funcDir = regexprep(funcDir(1).file,[funcDir(1).name '\.m'],'');
    annPath = fullfile(funcDir,'..','model');
end

% Column order matches the existing tsv headers exactly, so a derived
% file is not just equivalent but byte-for-byte comparable to a
% hand-maintained one.
writeAnnotationTsv(fullfile(annPath,'reactions.tsv'), model, 'rxns', 'rxnNames', ...
    {'bigg.reaction', 'ec-code', 'kegg.pathway', 'kegg.reaction', 'metanetx.reaction'}, ...
    {'rxnMiriams',    'eccodes', 'rxnMiriams',   'rxnMiriams',    'rxnMiriams'});

writeAnnotationTsv(fullfile(annPath,'metabolites.tsv'), model, 'mets', 'metNames', ...
    {'bigg.metabolite', 'chebi',      'kegg.compound', 'metanetx.chemical', 'smiles'}, ...
    {'metMiriams',      'metMiriams', 'metMiriams',     'metMiriams',       'metSmiles'});

writeGeneTsv(fullfile(annPath,'genes.tsv'), model);
end

function writeAnnotationTsv(filename, model, idField, nameField, columns, sourceFields)
% Write one id/name/<columns...> tsv. Each column reads from either a
% MIRIAM struct array field (matched by column name) or a plain
% per-entity scalar field (matched by having no MIRIAM entries under that
% name) -- sourceFields{c} names which model field backs columns{c}.
ids   = model.(idField);
names = model.(nameField);

fid = fopen(filename,'w');
closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'id\tname\t%s\n', strjoin(columns,'\t'));

for i = 1:numel(ids)
    cellValues = repmat({''}, 1, numel(columns));
    for c = 1:numel(columns)
        field = sourceFields{c};
        if ~isfield(model,field)
            continue
        end
        entry = model.(field){i};
        if isstruct(entry)
            % MIRIAM struct array: pick out this column's namespace.
            cellValues{c} = strjoin(entry.value(strcmp(entry.name, columns{c})), ';');
        elseif ~isempty(entry)
            % Plain scalar field (eccodes, metSmiles): the whole value.
            cellValues{c} = entry;
        end
    end
    fprintf(fid,'%s\t%s\t%s\n', ids{i}, names{i}, strjoin(cellValues,'\t'));
end
end

function writeGeneTsv(filename, model)
fid = fopen(filename,'w');
closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'id\tuniprot\n');

for i = 1:numel(model.genes)
    value = '';
    if isfield(model,'geneMiriams') && isstruct(model.geneMiriams{i})
        miriam = model.geneMiriams{i};
        value = strjoin(miriam.value(strcmp(miriam.name,'uniprot')), ';');
    end
    fprintf(fid,'%s\t%s\n', model.genes{i}, value);
end
end

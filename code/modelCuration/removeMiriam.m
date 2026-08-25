function model = removeMiriam(model,type,object,miriamName)
% removeMiriam
%   Removes all MIRIAM annotations of a given namespace from the specified
%   metabolites, reactions or genes. Complements RAVEN's editMiriam, which
%   can only add, fill or replace annotations, but not delete them.
%
%   model       model structure.
%   type        'met', 'rxn' or 'gene', dependent on which objects the
%               annotations should be removed from.
%   object      either a cell array of identifiers, a logical vector with
%               the same number of elements as the type (see above) in the
%               model, or a vector of indexes.
%   miriamName  string specifying the namespace of the identifiers that
%               should be removed, for instance 'kegg.compound'.
%
% Usage: model = removeMiriam(model,type,object,miriamName)

type = char(type);
if ~any(strcmp(type,{'met','rxn','gene'}))
    error('Invalid ''type'', should be ''met'', ''rxn'' or ''gene''.')
end
miriamName = char(miriamName);

if islogical(object)
    idx = transpose(find(object));
elseif isnumeric(object)
    idx = transpose(object(:));
else
    idx = transpose(getIndexes(model,object,[type 's']));
end

field = [type 'Miriams'];
if ~isfield(model,field)
    return
end

for i = idx
    miriam = model.(field){i};
    if isempty(miriam)
        continue
    end
    toKeep = ~strcmp(miriam.name,miriamName);
    if all(toKeep)
        continue
    end
    miriam.name = miriam.name(toKeep);
    miriam.value = miriam.value(toKeep);
    if isempty(miriam.name)
        model.(field){i} = [];
    else
        model.(field){i} = miriam;
    end
end
end

function model = saveYeastYaml(model,allowNoGrowth,deriveTsvs)
% saveYeastYaml
%   Save a curated model back to model/yeast-GEM.yml. Applies the
%   minimal_Y6 canonical medium, adds SBO terms, and checks aerobic and
%   anaerobic growth, then writes model/yeast-GEM.yml, with the
%   cross-reference annotation that lives in model/{reactions,
%   metabolites,genes}.tsv (yeast-GEM#379) and any deltaG fields
%   stripped back out -- deltaG is an estimated, not curator-verified
%   value, so it never goes into model/yeast-GEM.yml; call saveDeltaG
%   explicitly if you want to persist it to its own tsv files.
%
%   This is the function every curation script should call after editing
%   a model -- for producing the .xml/.txt/.xlsx/.mat files instead, see
%   commitYeastModel.
%
% Inputs:
%   model           (struct) model to save. Preferably RAVEN format,
%                   although COBRA format is also allowed, but some
%                   fields might be lost in the conversion.
%   allowNoGrowth   (bool, opt) if saving should be allowed whenever the
%                   model cannot grow, returning a warning (default
%                   true), otherwise will error.
%   deriveTsvs      (bool, opt) if true, model/reactions.tsv,
%                   model/metabolites.tsv and model/genes.tsv are
%                   overwritten from the model's current annotation via
%                   deriveAnnotationTsvs, so that an annotation edit made
%                   programmatically (not by hand-editing a tsv cell) is
%                   not silently dropped. Default false: a curator's
%                   hand-edited tsv is never overwritten unless this is
%                   explicitly requested.
%
% Output:
%   model   the model, with minimal_Y6 and SBO terms applied.
%
% Usage: model = saveYeastYaml(model,allowNoGrowth,deriveTsvs)

if nargin < 2 || isempty(allowNoGrowth)
    allowNoGrowth = true;
end
if nargin < 3 || isempty(deriveTsvs)
    deriveTsvs = false;
end
if ~(exist('ravenCobraWrapper.m','file')==2)
    error(['RAVEN cannot be found. See README.md for installation '...
        'instructions. RAVEN is required to make sure that the model '...
        'is stored in the correct file format for use in the '...
        'yeast-GEM GitHub repository'])
end

% Export as RAVEN format
if isfield(model,'rules')
    model = ravenCobraWrapper(model);
end

%Get and change to the script folder, as all folders are relative to this
%folder
scriptFolder = fileparts(which(mfilename));
currentDir = cd(scriptFolder);

%Set minimal media
cd modelCuration
model = minimal_Y6(model);
cd ..

%Update SBO terms in model:
cd missingFields
model = addSBOterms(model);
cd ..

%Check if model can grow:
checkGrowth(model,'aerobic',allowNoGrowth)
checkGrowth(model,'anaerobic',allowNoGrowth)

%Optionally re-derive the annotation tsvs from the model's own state,
%before stripping that same annotation back out of the .yml below.
if deriveTsvs
    deriveAnnotationTsvs(model,fullfile('..','model'));
end

%Update .yml model. Reaction, metabolite and gene cross-reference
%annotation (KEGG, BiGG, ChEBI, MetaNetX, EC codes, UniProt) lives in
%model/{reactions,metabolites,genes}.tsv instead (yeast-GEM#379), and
%deltaG lives in its own tsv files, so neither is re-embedded here.
leanModel = stripTsvAnnotation(model);
exportForGit(leanModel,'yeast-GEM','../model',{'yml'},false,false);

%Switch back to original folder
cd(currentDir)
end

%%
function model = stripTsvAnnotation(model)
%Remove the six cross-reference fields that live in
%model/{reactions,metabolites,genes}.tsv, and the deltaG fields that live
%in their own tsv files, from a copy of the model, so none of them are
%re-embedded in model/yeast-GEM.yml on every save (yeast-GEM#379). sbo
%and every other field is left untouched.
if isfield(model,'metDeltaG')
    model = rmfield(model,'metDeltaG');
end
if isfield(model,'rxnDeltaG')
    model = rmfield(model,'rxnDeltaG');
end
if isfield(model,'rxnMiriams')
    model.rxnMiriams = stripMiriamNames(model.rxnMiriams, ...
        {'bigg.reaction','kegg.pathway','kegg.reaction','metanetx.reaction'});
end
if isfield(model,'eccodes')
    model.eccodes = repmat({''}, size(model.eccodes));
end

if isfield(model,'metMiriams')
    model.metMiriams = stripMiriamNames(model.metMiriams, ...
        {'bigg.metabolite','chebi','kegg.compound','metanetx.chemical'});
end
if isfield(model,'metSmiles')
    model.metSmiles = repmat({''}, size(model.metSmiles));
end

if isfield(model,'geneMiriams')
    model.geneMiriams = stripMiriamNames(model.geneMiriams, {'uniprot'});
end
end

function miriams = stripMiriamNames(miriams, names)
%Remove any (name,value) pair whose name is in NAMES from every entity's
%MIRIAM struct, leaving other namespaces (e.g. sbo) untouched. An entity
%left with no annotation at all collapses to '', matching
%readYAMLmodel's own empty-entry convention.
for i = 1:numel(miriams)
    if ~isstruct(miriams{i})
        continue
    end
    keep = ~ismember(miriams{i}.name, names);
    if all(keep)
        continue
    elseif any(keep)
        miriams{i} = struct('name',{miriams{i}.name(keep)},'value',{miriams{i}.value(keep)});
    else
        miriams{i} = '';
    end
end
end

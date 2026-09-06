function model = loadYeastYaml(filename)
% loadYeastYaml
%   Load model/yeast-GEM.yml for curation, requires the RAVEN Toolbox.
%   Merges in the reaction, metabolite and gene cross-reference
%   annotation from model/reactions.tsv, model/metabolites.tsv and
%   model/genes.tsv (yeast-GEM#379).
%
%   Never carries deltaG: those are estimated, not curator-verified
%   values, so model/yeast-GEM.yml does not carry them and this strips
%   any that readYAMLmodel happens to find (e.g. from an older .yml).
%   Call loadDeltaG explicitly if you want them.
%
%   Loading model/yeast-GEM.xml does not need a yeast-GEM wrapper: a
%   generic loader (importModel, readCbModel or cobra.io.read_sbml_model)
%   already returns the complete model.
%
% Input:
%   filename    by default, the model is loaded from its location at
%               yeast-GEM/model/yeast-GEM.yml. An alternative .yml file
%               can be loaded if provided here (opt, default empty).
%
% Output:
%   model       the yeast-GEM model structure, with reaction, metabolite
%               and gene cross-reference annotation merged back in.
%
%   Usage: model = loadYeastYaml(filename)

funcDir = dbstack('-completenames');
funcDir = regexprep(funcDir(1).file,[funcDir(1).name '\.m'],'');

if nargin<1 || isempty(filename)
    filename = fullfile(funcDir,'..','model','yeast-GEM.yml');
end

model = readYAMLmodel(filename);

% Never surface deltaG from a load: it is an estimated, not
% curator-verified field (call loadDeltaG explicitly if you want it).
if isfield(model,'metDeltaG')
    model = rmfield(model,'metDeltaG');
end
if isfield(model,'rxnDeltaG')
    model = rmfield(model,'rxnDeltaG');
end

% The tsvs merged in here always come from this checkout's model/ folder,
% regardless of which .yml was loaded -- e.g. getEarlierModelVersion.m
% loads a git-shown revision into a temp filename, and still wants it
% annotated the same way a fresh default-path load would be.
annPath = fullfile(funcDir,'..','model');
if isfile(fullfile(annPath,'reactions.tsv'))
    model = annotateGEM(model, annPath);
end
end

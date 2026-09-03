function model = loadYeastYaml(filename)
% loadYeastYaml
%   Load model/yeast-GEM.yml for curation, requires the RAVEN Toolbox.
%   Merges in the reaction, metabolite and gene cross-reference
%   annotation from model/reactions.tsv, model/metabolites.tsv and
%   model/genes.tsv (yeast-GEM#379), and restores the deltaG fields.
%
%   Loading model/yeast-GEM.xml does not need a yeast-GEM wrapper: once
%   RAVEN/raven-toolbox support deltaG in SBML, a generic loader
%   (importModel, readCbModel or cobra.io.read_sbml_model) already
%   returns the complete model.
%
% Input:
%   filename    by default, the model is loaded from its location at
%               yeast-GEM/model/yeast-GEM.yml. An alternative .yml file
%               can be loaded if provided here (opt, default empty).
%
% Output:
%   model       the yeast-GEM model structure, with reaction, metabolite
%               and gene cross-reference annotation and deltaG fields
%               merged back in.
%
%   Usage: model = loadYeastYaml(filename)

funcDir = dbstack('-completenames');
funcDir = regexprep(funcDir(1).file,[funcDir(1).name '\.m'],'');

if nargin<1 || isempty(filename)
    filename = fullfile(funcDir,'..','model','yeast-GEM.yml');
end

model = readYAMLmodel(filename);

% The tsvs merged in here always come from this checkout's model/ folder,
% regardless of which .yml was loaded -- e.g. getEarlierModelVersion.m
% loads a git-shown revision into a temp filename, and still wants it
% annotated the same way a fresh default-path load would be.
annPath = fullfile(funcDir,'..','model');
if isfile(fullfile(annPath,'reactions.tsv'))
    model = annotateGEM(model, annPath);
end

currentDir = pwd;
cd(fullfile(funcDir,'missingFields'));
model = loadDeltaG(model);
cd(currentDir)
end

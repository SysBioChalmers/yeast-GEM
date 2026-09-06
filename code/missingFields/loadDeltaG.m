function model = loadDeltaG(model)
% loadDeltaG
%   Add metDeltaG and rxnDeltaG fields to a model, based on datafiles saved at
%   /data/databases (model_rxnDeltaG.tsv and model_metDeltaG.tsv). Metabolites
%   and reactions are matched by their identifiers (i.e. model.mets and
%   model.rxns). If changes are made that affect the identifiers or what
%   metabolites or reactions they refer to, the deltaG values will not be
%   correct.
%
%   These are estimated values, not curator-verified, so they are never
%   part of model/yeast-GEM.yml or the exported .xml/.txt/.xlsx/.mat --
%   call this explicitly if you want them (loadYeastYaml/saveYeastYaml
%   and commitYeastModel never call it).
%
% Input:
%   model   yeast-GEM without deltaG fields
%
% Output:
%   model   yeast-GEM with metDeltaG and rxnDeltaG fields
%
% Usage: model = loadDeltaG(model)

if isfield(model,'metDeltaG')
    disp('Existing metDeltaG field will be overwritten.')
else
    model.metDeltaG = nan(numel(model.mets),1);
end
if isfield(model,'rxnDeltaG')
    disp('Existing rxnDeltaG field will be overwritten.')
else
    model.rxnDeltaG = nan(numel(model.rxns),1);
end

metG = readtable('../../data/databases/model_metDeltaG.tsv', ...
    'FileType','text', 'Delimiter','\t');
rxnG = readtable('../../data/databases/model_rxnDeltaG.tsv', ...
    'FileType','text', 'Delimiter','\t');

[a,b] = ismember(model.mets,metG.id);
model.metDeltaG(a) = metG.deltaG(b(a));
if any(~a)
    fprintf(['Not all metabolite identifiers are matched to model_metDeltaG.tsv, the latter \n' ...
             'file might have to be supplemented with deltaG values for new metabolites.\n'])
end

[a,b] = ismember(model.rxns,rxnG.id);
model.rxnDeltaG(a) = rxnG.deltaG(b(a));
if any(~a)
    fprintf(['Not all reaction identifiers are matched to model_rxnDeltaG.tsv, the latter \n' ...
             'file might have to be supplemented with deltaG values for new reaction.\n'])
end
end

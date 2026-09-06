function model = saveDeltaG(model,verbose)
% saveDeltaG
%   Saves the metDeltaG and rxnDeltaG fields as tables to /data/databases/...
%   model_rxnDeltaG.tsv and /data/databases/model_metDeltaG.tsv. Call
%   loadDeltaG to reconstruct the metDeltaG and rxnDeltaG fields from
%   these files.
%
%   These are estimated values, not curator-verified, so they are never
%   part of model/yeast-GEM.yml or the exported .xml/.txt/.xlsx/.mat --
%   call this explicitly if you want to persist them (loadYeastYaml/
%   saveYeastYaml and commitYeastModel never call it).
%
% Input:
%   model       yeast-GEM with deltaG fields
%   verbose     true or false
%
% Output:
%   model   yeast-GEM with deltaG fields
%
% Usage: model = saveDeltaG(model,verbose)

if nargin<2
    verbose=false;
end
if ~isfield(model,'metDeltaG')
    if verbose
        disp('No metDeltaG field found, model_metDeltaG.tsv will not be changed.')
    end
else
    metG = table(model.mets, model.metDeltaG, 'VariableNames',{'id','deltaG'});
    writetable(metG,'../../data/databases/model_metDeltaG.tsv', ...
        'FileType','text', 'Delimiter','\t');
end
if ~isfield(model,'rxnDeltaG')
    if verbose
        disp('No rxnDeltaG field found, model_rxnDeltaG.tsv will not be changed')
    end
else
    rxnG = table(model.rxns, model.rxnDeltaG, 'VariableNames',{'id','deltaG'});
    writetable(rxnG,'../../data/databases/model_rxnDeltaG.tsv', ...
        'FileType','text', 'Delimiter','\t');
end
end

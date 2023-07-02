function model = saveDeltaG(model,verbose)
% saveDeltaG
%   Saves the metDeltaG and rxnDeltaG fields as tables to /data/databases/...
%   model_rxnDeltaG.csv and /data/databases/model_metDeltaG.csv. When
%   loadYeastModel is run, these tables will be read to reconstruct the
%   metDeltaG and rxnDeltaG fields.
%
% Input:
%   model       yeast-GEM with deltaG fields
%   verbose     true or false
%
% Output:
%   model   yeast-GEM with metDeltaG and rxnDeltaG fields
%
% Usage: model = saveDeltaG(model,verbose)

if nargin<2
    verbose=false;
end
if ~isfield(model,'metDeltaG')
    if verbose
        disp('No metDeltaG field found, model_metDeltaG.csv will not be changed.')
    end
else
    metG = array2table([model.mets, num2cell(model.metDeltaG)]);
    writetable(metG,'../../data/databases/model_metDeltaG.csv');
end
if ~isfield(model,'rxnDeltaG')
    if verbose
        disp('No rxnDeltaG field found, model_rxnDeltaG.csv will not be changed')
    end
else
    rxnG = array2table([model.rxns, num2cell(model.rxnDeltaG)]);
    writetable(rxnG,'../../data/databases/model_rxnDeltaG.csv');
end
end

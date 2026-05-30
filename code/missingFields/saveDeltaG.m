function model = saveDeltaG(model, verbose)
% saveDeltaG  yeast-GEM shim — delegates to RAVEN's saveDeltaGtoCSV.
%
%   Persists model.metDeltaG and model.rxnDeltaG to the project
%   CSVs at data/databases/model_metDeltaG.csv and
%   data/databases/model_rxnDeltaG.csv. Returns the model unchanged
%   (kept as an output for backward compatibility with callers that
%   chain saveDeltaG into a pipeline).
%
% Usage: model = saveDeltaG(model)
%        model = saveDeltaG(model, verbose)

if nargin < 2
    verbose = false;
end

funcDir = fileparts(mfilename('fullpath'));
metCsv  = fullfile(funcDir, '..', '..', 'data', 'databases', 'model_metDeltaG.csv');
rxnCsv  = fullfile(funcDir, '..', '..', 'data', 'databases', 'model_rxnDeltaG.csv');
saveDeltaGtoCSV(model, metCsv, rxnCsv, verbose);
end

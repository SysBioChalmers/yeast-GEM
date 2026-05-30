function model = loadDeltaG(model)
% loadDeltaG  yeast-GEM shim — delegates to RAVEN's loadDeltaGfromCSV.
%
%   Populates model.metDeltaG and model.rxnDeltaG from the project
%   CSVs at data/databases/model_metDeltaG.csv and
%   data/databases/model_rxnDeltaG.csv. Paths are resolved relative to
%   this file so the function works from any cwd.
%
% Usage: model = loadDeltaG(model)

funcDir = fileparts(mfilename('fullpath'));
metCsv  = fullfile(funcDir, '..', '..', 'data', 'databases', 'model_metDeltaG.csv');
rxnCsv  = fullfile(funcDir, '..', '..', 'data', 'databases', 'model_rxnDeltaG.csv');
model   = loadDeltaGfromCSV(model, metCsv, rxnCsv);
end

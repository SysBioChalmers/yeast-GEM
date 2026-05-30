function runPhase6Curation(yeastGemPath, outXml)
% runPhase6Curation  Apply the v8_6_3 VolPolyP curation TSVs and export.
%
%   Used for the Python-vs-MATLAB curation-engine parity check.

warning('off','all');
restoredefaultpath; rehash toolboxcache;
addpath('/opt/gurobi1301/linux64/matlab');
addpath(genpath('/home/eduardk/github/RAVEN'));
addpath(fullfile(yeastGemPath, 'code'));
addpath(fullfile(yeastGemPath, 'code', 'modelCuration'));
addpath(fullfile(yeastGemPath, 'code', 'otherChanges'));
addpath(fullfile(yeastGemPath, 'code', 'missingFields'));

dataDir = fullfile(yeastGemPath, 'data', 'modelCuration', 'v8_7_0');
mets   = fullfile(dataDir, 'DBnewRxnsMets.tsv');
genes  = fullfile(dataDir, 'DBnewRxnsGenes.tsv');
rxns   = fullfile(dataDir, 'DBnewRxnsRxns.tsv');
coeffs = fullfile(dataDir, 'DBnewRxnsCoeffs.tsv');

model = loadYeastModel;
model = curateMetsRxnsGenes(model, mets, genes, coeffs, rxns);
exportModel(model, outXml);
fprintf('Wrote %s\n', outXml);
end

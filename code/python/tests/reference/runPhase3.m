function runPhase3(yeastGemPath, outFile, funcName)
% runPhase3  Drive saveYeastModel / commitYeastModel for phase-3 verification.
%
%   model = loadYeastModel
%   model = <funcName>(model, false, true, false)  % no README update,
%                                                  % allow no growth, no binary
%   copyfile model/yeast-GEM.xml -> outFile

warning('off','all');
restoredefaultpath; rehash toolboxcache;
addpath(genpath('/home/eduardk/github/RAVEN'));
addpath(fullfile(yeastGemPath, 'code'));
addpath(fullfile(yeastGemPath, 'code', 'modelCuration'));
addpath(fullfile(yeastGemPath, 'code', 'otherChanges'));
addpath(fullfile(yeastGemPath, 'code', 'missingFields'));

model = loadYeastModel;
fcn = str2func(funcName);
fcn(model, false, true, false);

src = fullfile(yeastGemPath, 'model', 'yeast-GEM.xml');
copyfile(src, outFile);
fprintf('Wrote %s\n', outFile);
end

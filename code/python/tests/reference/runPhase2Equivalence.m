function runPhase2Equivalence(yeastGemPath, outDir)
% runPhase2Equivalence  Apply the four phase-2 conditions and save outputs.
%
%   Drives the MATLAB-side equivalence check for the phase-2 refactor.
%   For each of the four conditions:
%     - Saves the full RAVEN model struct as <name>.mat (always works,
%       even for infeasible bound states).
%     - Additionally exports SBML to <name>.xml when the bounds are
%       feasible (lb <= ub everywhere); skipped otherwise with a note.
%
%   The .mat files are compared via comparePhase2.m; the .xml files are
%   compared via `python -m yeastgem.compare` (level-1 semantic gate).

warning('off','all');
restoredefaultpath; rehash toolboxcache;
addpath(genpath('/home/eduardk/github/RAVEN'));
addpath(fullfile(yeastGemPath, 'code'));
addpath(fullfile(yeastGemPath, 'code', 'modelCuration'));
addpath(fullfile(yeastGemPath, 'code', 'otherChanges'));
addpath(fullfile(yeastGemPath, 'code', 'missingFields'));

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

conds = {'minimal_Y6', 'anaerobicModel', 'glycineNitrogenSource', 'nitrogenLimitation'};

for i = 1:numel(conds)
    name = conds{i};
    fprintf('=== %s ===\n', name);
    model = loadYeastModel;
    fcn = str2func(name);
    model = fcn(model); %#ok<NASGU>

    matFile = fullfile(outDir, [name '.mat']);
    save(matFile, 'model', '-v7');
    fprintf('Wrote %s\n', matFile);

    feasible = all(model.lb <= model.ub);
    if feasible
        xmlFile = fullfile(outDir, [name '.xml']);
        exportModel(model, xmlFile);
        fprintf('Wrote %s\n', xmlFile);
    else
        nBad = sum(model.lb > model.ub);
        fprintf('SKIP %s.xml (%d reactions have lb > ub; SBML export would fail)\n', ...
            name, nBad);
    end
end
fprintf('=== All conditions applied successfully ===\n');
end

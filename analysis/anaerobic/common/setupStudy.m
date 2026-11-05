function paths = setupStudy(ravenVersion)
% setupStudy  Clean MATLAB path for the anaerobic analyses.
%   paths = setupStudy('3.0.0-rc1')
%
%   Puts one RAVEN checkout, the Gurobi MATLAB interface, the yeast-GEM code
%   and the analysis functions on a fresh path, sets RAVEN's solver to Gurobi,
%   and shadows savepath so the pathdef.m of the user is never rewritten.
%   ravenVersion names the RAVEN version that the results were made with; it is
%   only printed.
%
%   Environment variables:
%     RAVEN_ROOT       the RAVEN checkout (required)
%     GUROBI_HOME      the Gurobi installation; its MATLAB interface is in matlab/
%     YEAST_GEM_ROOT   the yeast-GEM checkout whose code, data and git tags are
%                      used (default: the repository that contains this folder)
%     ANAEROBIC_RESULTS  folder for the result files (default: results/ in this folder)
%     MATLAB_PREFDIR   use a folder of its own, so that RAVEN's setpref stays out
%                      of the user's MATLAB preferences

study = fileparts(fileparts(mfilename('fullpath')));
paths.study  = study;
paths.yeast  = yeastRoot();
paths.raven  = getenv('RAVEN_ROOT');
paths.gurobi = fullfile(getenv('GUROBI_HOME'), 'matlab');
paths.results = getenv('ANAEROBIC_RESULTS');
if isempty(paths.results), paths.results = fullfile(study, 'results'); end
if ~isfolder(paths.results), mkdir(paths.results); end

assert(isfolder(paths.raven), 'Set RAVEN_ROOT to the RAVEN checkout (RAVEN %s was used).', ravenVersion);
assert(isfolder(paths.gurobi), 'Gurobi MATLAB interface not found; set GUROBI_HOME');
if isempty(getenv('MATLAB_PREFDIR'))
    warning('setupStudy:prefdir', 'MATLAB_PREFDIR is not set; RAVEN will write to your MATLAB preferences.');
end

restoredefaultpath;
ravenPath = strsplit(genpath(paths.raven), pathsep);
ravenPath = ravenPath(~contains(ravenPath, [filesep '.git']) & ~cellfun(@isempty, ravenPath));
addpath(strjoin(ravenPath, pathsep));
addpath(paths.gurobi);
addpath(fullfile(paths.yeast, 'code'));
addpath(fullfile(paths.yeast, 'code', 'otherChanges'));
folders = strsplit(genpath(study), pathsep);
addpath(strjoin(folders(~contains(folders, [filesep 'results']) & ~cellfun(@isempty, folders)), pathsep));
addpath(fullfile(study, 'common', 'guard'));   % last, so it is first on the path

warning('off', 'yeastGEM:anaerobicModel:deprecated');
setRavenSolver('gurobi');
fprintf('RAVEN (%s) from %s\n', ravenVersion, which('solveLP'));
fprintf('Gurobi from %s\n', which('gurobi'));
end

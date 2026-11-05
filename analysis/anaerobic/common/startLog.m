function cleanup = startLog(paths, name)
% startLog  Copy everything printed to the command window to results/<name>.txt
%   until the returned object is cleared (end of the calling function).
if ~isfolder(paths.results), mkdir(paths.results); end
file = fullfile(paths.results, [name '.txt']);
if exist(file, 'file'), delete(file); end
diary(file);
cleanup = onCleanup(@() diary('off'));
end

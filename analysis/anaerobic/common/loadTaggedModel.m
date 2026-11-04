function model = loadTaggedModel(yeastDir, version)
% loadTaggedModel  Load model/yeast-GEM.yml from a yeast-GEM release tag.
%   Same as getEarlierModelVersion, but writes the YAML to a file private to
%   this MATLAB process, so jobs running at the same time in the same folder
%   cannot read each other's model.

tmp = [tempname() '.yml'];
cleanup = onCleanup(@() delete(tmp));
cmd = sprintf('git -C "%s" show refs/tags/v%s:model/yeast-GEM.yml > "%s"', yeastDir, version, tmp);
status = system(cmd);
assert(status == 0, 'git show failed for version %s', version);
model = loadYeastModel(tmp);
loaded = regexprep(model.id, 'yeastGEM_v?', '');
assert(strcmp(loaded, version), 'Loaded model %s, expected %s', loaded, version);
fprintf('Loaded model version: ''%s''\n', loaded);
end

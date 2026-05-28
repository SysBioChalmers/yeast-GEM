function out = readYAML(filename)
% readYAML  Read an arbitrary YAML file into a MATLAB struct.
%
%   Used by applyCondition / applyIDs to consume the data files under
%   data/conditions/ and data/yeastgem/ that drive the data-as-code
%   refactor (PORTING_PLAN.md, phase 2).
%
%   Implementation: delegates to Python's `yaml.safe_load`, then
%   recursively converts the returned py.dict / py.list to native
%   MATLAB struct / cell. Requires a working MATLAB-Python bridge
%   (i.e. `pyversion` set) with `pyyaml` installed:
%
%       pip install pyyaml         % from the MATLAB-linked Python env
%
%   This avoids both writing a MATLAB YAML parser and depending on
%   `yamlread` (MATLAB R2024a+). RAVEN's own YAML reader handles only
%   the cobra-format model layout, not arbitrary YAML.
%
% Usage: cfg = readYAML('data/conditions/anaerobic.yml')

if ~isfile(filename)
    error('readYAML:fileNotFound', 'File not found: %s', filename);
end

try
    py.importlib.import_module('yaml');
catch ME
    error('readYAML:pyyamlMissing', ...
        ['pyyaml is required to read the data files. Install it in ' ...
         'your MATLAB-linked Python environment (`pip install pyyaml`).' ...
         '\nUnderlying error: %s'], ME.message);
end

f = py.builtins.open(filename, 'r');
try
    data = py.yaml.safe_load(f);
finally
    f.close();
end

out = pyToMatlab(data);
end


function v = pyToMatlab(obj)
% Recursively convert Python objects from pyyaml into MATLAB types.
if isa(obj, 'py.NoneType')
    v = [];
elseif isa(obj, 'py.bool')
    v = logical(obj);
elseif isa(obj, 'py.int') || isa(obj, 'py.float')
    v = double(obj);
elseif isa(obj, 'py.str')
    v = char(obj);
elseif isa(obj, 'py.dict')
    v = struct();
    keys = cell(py.list(obj.keys()));
    vals = cell(py.list(obj.values()));
    for i = 1:numel(keys)
        key = char(keys{i});
        v.(matlabFieldName(key)) = pyToMatlab(vals{i});
    end
elseif isa(obj, 'py.list') || isa(obj, 'py.tuple')
    cells = cell(obj);
    v = cell(numel(cells), 1);
    for i = 1:numel(cells)
        v{i} = pyToMatlab(cells{i});
    end
else
    % Fallback: best-effort scalar
    v = obj;
end
end


function name = matlabFieldName(key)
% Sanitise a YAML key into a valid MATLAB field name.
% Replaces non-alphanumeric characters with underscores; prefixes a
% digit-starting key with `f_`.
name = regexprep(key, '[^A-Za-z0-9_]', '_');
if isempty(name) || ~isstrprop(name(1), 'alpha')
    name = ['f_' name];
end
end

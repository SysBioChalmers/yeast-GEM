function saveYeastModel(model,upDATE,allowNoGrowth,binaryFiles)
% saveYeastModel  DEPRECATED — renamed to commitYeastModel.
%
%   `saveYeastModel` implied a casual write, but the function is the
%   heavy release pipeline you run before opening a curation PR. It is
%   renamed to `commitYeastModel` to make that workflow explicit. The
%   docstring of `commitYeastModel` clarifies that it does NOT perform
%   `git commit`; it prepares the artifacts so the next `git commit`
%   captures a coherent release-ready state.
%
%   This shim forwards all arguments to `commitYeastModel` and emits a
%   deprecation warning. It will be removed at the next minor version
%   bump after the rename ships.
%
% Usage: saveYeastModel(model,upDATE,allowNoGrowth,binaryFiles)

warning('yeastGEM:saveYeastModelDeprecated', ...
    ['saveYeastModel is deprecated; use commitYeastModel instead. ' ...
     'See code/python/PORTING_PLAN.md (phase 3) for the rename rationale.']);

if nargin < 2
    commitYeastModel(model);
elseif nargin < 3
    commitYeastModel(model, upDATE);
elseif nargin < 4
    commitYeastModel(model, upDATE, allowNoGrowth);
else
    commitYeastModel(model, upDATE, allowNoGrowth, binaryFiles);
end
end

function model = saveYeastModel(model,upDATE,allowNoGrowth,binaryFiles)
% saveYeastModel
%   DEPRECATED. Use saveYeastYaml for a routine curation save (writes
%   model/yeast-GEM.yml only), or commitYeastModel for a local binary
%   build or release (writes .xml/.txt/.xlsx/.mat).
%
%   Kept as a shim for external code that still calls this name. Forwards
%   to commitYeastModel, the closest behavioural match: unlike
%   saveYeastYaml, the original saveYeastModel always wrote .xml/.txt (and
%   .xlsx/.mat when binaryFiles was true), never just the .yml.
%
% Inputs:
%   model           (struct) model to save. Preferably RAVEN format,
%                   although COBRA format is also allowed, but some
%                   fields might be lost in the conversion.
%   upDATE          unused (kept for call-signature compatibility).
%   allowNoGrowth   (bool, opt) if saving should be allowed whenever the
%                   model cannot grow, returning a warning (default true),
%                   otherwise will error.
%   binaryFiles     (bool, opt) if the model should also be stored in
%                   binary file formats (= xlsx and mat).
%
% Usage: model = saveYeastModel(model,upDATE,allowNoGrowth,binaryFiles)

warning(['saveYeastModel is deprecated; use saveYeastYaml for a routine '...
    'curation save, or commitYeastModel to write .xml/.txt/.xlsx/.mat.'])

if nargin < 3 || isempty(allowNoGrowth)
    allowNoGrowth = true;
end
if nargin < 4 || isempty(binaryFiles)
    binaryFiles = false;
end

if binaryFiles
    formats = {'xml','txt','xlsx','mat'};
else
    formats = {'xml','txt'};
end

model = commitYeastModel(model,formats,allowNoGrowth);
end

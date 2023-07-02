function model = loadYeastModel
% loadYeastModel
%   Load the yeast-GEM in a MATLAB environment, requires the RAVEN Toolbox
%   and libSBML to be installed. By default, the model is provided in RAVEN
%   format by default, but can alternatively be converted to COBRA format.
%   If RAVEN is not installed, but COBRA is, then the model will instead
%   load the model via COBRA after throwing a warning that RAVEN is
%   recommended.
%
%   Input:
%       cobra       true if the model should be provided as COBRA format,
%                   false for RAVEN format (optional, default false)
%
%   Output:
%       model       the yeast-GEM model structure
%
%   Usage: model = loadYeastModel(cobra)

scriptFolder = fileparts(which(mfilename));
currentDir = cd(scriptFolder);
cd(currentDir)
if ~(exist('ravenCobraWrapper.m','file')==2)
    if exist('readCbModel.m','file')==2
        warning(['RAVEN cannot be found. yeast-GEM will instead be loaded in '...
            'COBRA format.\n\nNote that it is recommended to have RAVEN '...
            'installed, especially when curating yeast-GEM (see README.md for '...
            'more info).%s'],'')
        model = readCbModel('../model/yeast-GEM.xml');
    else
        error(['RAVEN cannot be found. See README.md for installation '...
            'instructions.'])
    end
else
    model = importModel('../model/yeast-GEM.xml');
    if nargin>0 && cobra==true
        model = ravenCobraWrapper(model);
    end
end
cd missingFields
model = loadDeltaG(model);
cd(currentDir)
end

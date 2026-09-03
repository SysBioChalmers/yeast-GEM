function checkGrowth(model,condition,allowNoGrowth)
% checkGrowth
%   Check if the model can grow under a given condition (aerobic or
%   anaerobic) using RAVEN. Will either return warnings or errors
%   depending on allowNoGrowth.
%
%   Extracted from saveYeastModel.m so that saveYeastYaml.m and
%   commitYeastModel.m can both run the same aerobic/anaerobic gate
%   without duplicating it.
%
% Input:
%   model           model structure to check.
%   condition       'aerobic' or 'anaerobic'.
%   allowNoGrowth   if true, a no-growth or infeasible result is a
%                   warning; if false, it is an error.
%
% Usage: checkGrowth(model,condition,allowNoGrowth)

if strcmp(condition,'anaerobic')
    cd otherChanges
    model = anaerobicModel(model);
    cd ..
end
try
    xPos = strcmp(model.rxnNames,'growth');
    sol  = solveLP(model);
    if sol.x(xPos) < 1e-6
        dispText = ['The model is not able to support growth under ' ...
            condition ' conditions. Please ensure the model can grow'];
    end
catch
    dispText = ['The model yields an infeasible simulation using RAVEN ' ...
        'under ' condition ' conditions. Please ensure the model ' ...
        'can be simulated with RAVEN'];
end

if exist('dispText','var')
    if allowNoGrowth
        warning([dispText ' before opening a PR.'])
    else
        error([dispText ' before committing.'])
    end
end
end

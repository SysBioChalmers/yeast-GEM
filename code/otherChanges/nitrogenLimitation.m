function model = nitrogenLimitation(model)
% nitrogenLimitation  Convert model to nitrogen-limiting conditions.
%
%   Now a thin shim around applyCondition('nitrogen_limitation'); the
%   bound changes live in data/conditions/nitrogen_limitation.yml. See
%   code/python/PORTING_PLAN.md (phase 2) for the data-as-code refactor.
%
%   Reference: doi:10.1128/EC.2.5.827-829.2003.
%
% Usage: model = nitrogenLimitation(model)

model = applyCondition(model, 'nitrogen_limitation');
end

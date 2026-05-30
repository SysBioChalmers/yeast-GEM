function model = glycineNitrogenSource(model)
% glycineNitrogenSource  Convert model to glycine-as-N-source.
%
%   Now a thin shim around applyCondition('glycine_nitrogen'); the bound
%   changes live in data/conditions/glycine_nitrogen.yml. See
%   code/python/PORTING_PLAN.md (phase 2) for the data-as-code refactor.
%
%   References:
%       doi:10.1111/j.1567-1364.2002.tb00069.x
%       doi:10.1074/jbc.274.15.10523
%       doi:10.1128/EC.2.5.827-829.2003
%
% Usage: model = glycineNitrogenSource(model)

model = applyYeastCondition(model, 'glycine_nitrogen');
end

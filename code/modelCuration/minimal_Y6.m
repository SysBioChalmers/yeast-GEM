function model = minimal_Y6(model)
% minimal_Y6  Set Y6 minimal-media exchange bounds.
%
%   Now a thin shim around applyCondition('minimal_Y6'); the bound
%   changes live in data/conditions/minimal_Y6.yml. See
%   code/python/PORTING_PLAN.md (phase 2) for the data-as-code refactor.
%
%   Reference: doi:10.1371/journal.pcbi.1004530.
%
% Usage: model = minimal_Y6(model)

model = applyYeastCondition(model, 'minimal_Y6');
end

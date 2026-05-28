function model = anaerobicModel(model)
% anaerobicModel  Constrain yeast-GEM to anaerobic conditions.
%
%   Now a thin shim around applyCondition('anaerobic'); the cofactor
%   pseudoreaction edits, amino-acid ratio switch, biomass
%   stoichiometry delta and bound changes live in
%   data/conditions/anaerobic.yml. See code/python/PORTING_PLAN.md
%   (phase 2) for the data-as-code refactor.
%
%   This function was last updated as part of release v9.1.0.
%
% Usage: model = anaerobicModel(model)

model = applyCondition(model, 'anaerobic');
end

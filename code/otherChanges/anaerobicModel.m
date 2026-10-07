function model = anaerobicModel(model)
% anaerobicModel
%   Deprecated. Applies the anaerobic constraints of yeast-GEM releases
%   before 9.1.0 (same as anaerobicModelOld) and warns. The anaerobic
%   constraints curated since yeast-GEM 9.1.0 are applied by
%   applyEnvironment(model,'anaerobic').
%
% Input:
%   model           yeast-GEM model structure, which is aerobic by default
%
% Output:
%   model           model structure with the pre-9.1.0 anaerobic constraints
%
% Usage: model = anaerobicModel(model)

warning('yeastGEM:anaerobicModel:deprecated', ...
    ['anaerobicModel is deprecated and applies the anaerobic constraints of ' ...
     'yeast-GEM releases before 9.1.0. Use applyEnvironment(model,''anaerobic'') ' ...
     'for the anaerobic constraints curated since yeast-GEM 9.1.0.']);
model = anaerobicModelOld(model);
end

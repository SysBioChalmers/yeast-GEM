function model = minimal_Y6(model)
% minimal_Y6
%   Minimal glucose medium (ammonium, glucose, oxygen, phosphate, sulphate
%   and trace elements), from doi:10.1371/journal.pcbi.1004530. Kept for
%   backwards compatibility: same as applyEnvironment(model,'minimal_Y6');
%   the constraints are in data/conditions/minimal_Y6.yml.
%
%   Usage: model = minimal_Y6(model)

model = applyEnvironment(model,'minimal_Y6');
end

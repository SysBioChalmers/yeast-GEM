function model = applyYeastCondition(model, name)
% applyYeastCondition  Apply a named yeast-GEM condition preset to the model.
%
%   Yeast-specific wrapper around RAVEN's generic applyCondition. This
%   function:
%       1. Resolves `name` to a YAML file under `data/conditions/`.
%       2. Applies the yeast-specific `amino_acid_ratio` step
%          (via changeAminoAcidRatio) when present in the YAML.
%       3. Hands the parsed condition to RAVEN's `applyCondition` for
%          the generic prelude / cofactor / biomass-delta / bounds /
%          uptake-count steps.
%
%   Available presets (data/conditions/<name>.yml):
%       'minimal_Y6'         minimal media (replaces minimal_Y6.m)
%       'anaerobic'          anaerobic conditions (replaces anaerobicModel.m)
%       'glycine_nitrogen'   glycine as sole N source
%       'nitrogen_limitation' N-limited
%
%   Requires RAVEN with the applyCondition / parseYAML helpers (commit
%   on the feat/yeast-gem-shared branch or any later release that
%   incorporates them).
%
% Usage: model = applyYeastCondition(model, 'anaerobic')

funcDir = fileparts(mfilename('fullpath'));
yamlPath = fullfile(funcDir, '..', 'data', 'conditions', [name '.yml']);
if ~isfile(yamlPath)
    error('applyYeastCondition:unknownCondition', ...
        'No such condition: %s (looked for %s)', name, yamlPath);
end
cond = parseYAML(yamlPath);

% Yeast-specific pre-step: amino_acid_ratio rewrites the protein
% pseudoreaction's stoichiometry from data/physiology/. The generic
% applyCondition silently ignores this field.
if isfield(cond, 'amino_acid_ratio')
    aerobic = strcmp(cond.amino_acid_ratio, 'aerobic');
    model = changeAminoAcidRatio(model, aerobic);
end

% Generic mechanism (provided by RAVEN).
model = applyCondition(model, cond);
end

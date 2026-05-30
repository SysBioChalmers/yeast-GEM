function ids = applyIDs()
% applyIDs  Load the canonical yeast-GEM identifiers from data/yeastgem/ids.yml.
%
%   ids = applyIDs() returns a struct with fields:
%       biomass_rxn          string
%       protein_rxn          string
%       cofactor_rxn         string
%       proton_met           string
%       pseudoreaction_names struct (component -> name)
%       gam_cofactors        cell array of strings
%
%   This is the data-driven replacement for the hardcoded IDs that
%   used to live in functions like changeGAM.m, rescalePseudoReaction.m,
%   sumBioMass.m. Those functions are kept as legacy shims; new code
%   should call applyIDs and read from the returned struct.
%
%   Requires RAVEN's parseYAML (any RAVEN release ≥ the commit that
%   added io/parseYAML.m, currently the feat/yeast-gem-shared branch).
%
% Usage: ids = applyIDs()

funcDir = fileparts(mfilename('fullpath'));
yamlPath = fullfile(funcDir, '..', 'data', 'yeastgem', 'ids.yml');
ids = parseYAML(yamlPath);
end

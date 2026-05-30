function cfg = yeastBiomassConfig()
% yeastBiomassConfig
%   Build the biomassConfig struct that RAVEN's biomass helpers expect,
%   populated from the canonical yeast-GEM IDs file
%   (data/yeastgem/ids.yml). All yeast-GEM-side biomass shims
%   (sumBioMass, scaleBioMass, rescalePseudoReaction, changeGAM) use
%   this so the IDs live in one place.
%
%   Output:
%       cfg     struct with the shape RAVEN's getBiomassFractions /
%               scaleBiomass*/ setGAM consume:
%                 .biomass_rxn          rxn id (string)
%                 .proton_met           met id (string)
%                 .components           cell array of structs:
%                                         .name
%                                         .pseudoreaction_name
%                                         .mass_strategy
%                 .gam_cofactors        cell array of met NAMES used
%                                       by changeGAM
%
% Usage: cfg = yeastBiomassConfig()

ids = applyIDs();

cfg.biomass_rxn = ids.biomass_rxn;
cfg.proton_met = ids.proton_met;
cfg.gam_cofactors = ids.gam_cofactors;

cfg.components = cell(numel(ids.biomass_components), 1);
for i = 1:numel(ids.biomass_components)
    comp = ids.biomass_components{i};
    cfg.components{i} = struct( ...
        'name', comp.name, ...
        'pseudoreaction_name', ids.pseudoreaction_names.(comp.name), ...
        'mass_strategy', comp.mass_strategy);
end
end

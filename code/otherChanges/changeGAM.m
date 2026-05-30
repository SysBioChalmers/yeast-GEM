function model = changeGAM(model, GAM, NGAM)
% changeGAM  yeast-GEM shim — delegates to RAVEN's setGAM.
%
%   Sets the GAM coefficient on the yeast-GEM biomass pseudoreaction
%   for the metabolites listed under `gam_cofactors` in
%   data/yeastgem/ids.yml (ATP, ADP, H2O, H+, phosphate by default).
%   If NGAM is supplied, the 'non-growth associated maintenance
%   reaction' is fixed to that flux.
%
% Usage: model = changeGAM(model, GAM)
%        model = changeGAM(model, GAM, NGAM)

cfg = yeastBiomassConfig();

if nargin > 2
    ngamPos = strcmp(model.rxnNames, 'non-growth associated maintenance reaction');
    ngamRxn = model.rxns{ngamPos};
    model = setGAM(model, GAM, cfg.biomass_rxn, cfg.gam_cofactors, ngamRxn, NGAM);
else
    model = setGAM(model, GAM, cfg.biomass_rxn, cfg.gam_cofactors);
end
end

function model = rescalePseudoReaction(model, metName, f)
% rescalePseudoReaction  yeast-GEM shim — delegates to RAVEN's
% scaleBiomassPseudoreaction, plus a yeast-only lipid aggregation.
%
%   metName 'lipid' rescales both 'lipid backbone' and 'lipid chain' —
%   yeast-GEM keeps lipid mass as backbone + chain, so users still
%   address it as a single component. 'lipid backbone' and 'lipid
%   chain' are handled directly here because the model.metNames for
%   their products contain a space and don't match the underscore-
%   compatible cfg.components{i}.name keys ('lipid_backbone', etc.)
%   that the RAVEN helper uses for product-side detection. Every
%   other component name is forwarded to RAVEN.
%
% Usage: model = rescalePseudoReaction(model, metName, f)

if strcmp(metName, 'lipid')
    model = rescalePseudoReaction(model, 'lipid backbone', f);
    model = rescalePseudoReaction(model, 'lipid chain', f);
    return;
end

cfg = yeastBiomassConfig();

if strcmp(metName, 'lipid backbone') || strcmp(metName, 'lipid chain')
    rxnName = [metName ' pseudoreaction'];
    rxnPos  = find(strcmp(model.rxnNames, rxnName));
    if isempty(rxnPos)
        return;
    end
    for i = 1:length(model.mets)
        S_ir   = model.S(i, rxnPos);
        isProd = strcmp(model.metNames{i}, metName);
        if S_ir ~= 0 && ~isProd
            model.S(i, rxnPos) = f * S_ir;
        end
    end
    Hc = find(strcmp(model.mets, cfg.proton_met));
    model.S(Hc, rxnPos) = 0;
    model.S(Hc, rxnPos) = -sum(model.S(:, rxnPos) .* model.metCharges, 'omitnan');
    return;
end

model = scaleBiomassPseudoreaction(model, cfg, metName, f);
end

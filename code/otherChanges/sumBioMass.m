function [X, P, C, R, D, L, I, F] = sumBioMass(model, dispOutput)
% sumBioMass  yeast-GEM shim — delegates to RAVEN's getBiomassFractions.
%
%   Calls getBiomassFractions(model, yeastBiomassConfig()) and unpacks
%   the resulting struct into the legacy (X, P, C, R, D, L, I, F)
%   outputs (total / protein / carbohydrate / RNA / DNA / lipid
%   backbone / ion / cofactor). yeast-GEM callers (scaleBioMass,
%   changeAminoAcidRatio, the v8_*/v9_* curation scripts) keep their
%   existing signatures.
%
%   Requires RAVEN ≥ the commit that added core/getBiomassFractions.m.
%
% Usage: [X, P, C, R, D, L, I, F] = sumBioMass(model, dispOutput)

if nargin < 2
    dispOutput = true;
end

cfg = yeastBiomassConfig();
fractions = getBiomassFractions(model, cfg);

P = fractions.protein;
C = fractions.carbohydrate;
R = fractions.RNA;
D = fractions.DNA;
L = fractions.lipid_backbone;
I = fractions.ion;
F = fractions.cofactor;
X = fractions.total;

if dispOutput
    disp(['P -> ' num2str(P) ' g/gDW'])
    disp(['C -> ' num2str(C) ' g/gDW'])
    disp(['R -> ' num2str(R) ' g/gDW'])
    disp(['D -> ' num2str(D) ' g/gDW'])
    disp(['L -> ' num2str(L) ' g/gDW'])
    disp(['I -> ' num2str(I) ' g/gDW'])
    disp(['F -> ' num2str(F) ' g/gDW'])
    disp(['X -> ' num2str(X) ' gDW/gDW'])
    sol = solveLP(model, 1);
    disp(['Growth = ' num2str(sol.f) ' 1/h'])
    disp(' ')
end
end

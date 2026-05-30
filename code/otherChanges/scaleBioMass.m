function model = scaleBioMass(model, component, new_value, balance_out, dispOutput)
% scaleBioMass  yeast-GEM shim — delegates to RAVEN's scaleBiomassFraction.
%
%   Scales `component` to `new_value` g/gDW, optionally adjusting
%   `balance_out` so the biomass total stays at 1 g/gDW. yeast-GEM's
%   'lipid' aggregation (rescale both backbone + chain in lock-step)
%   is preserved via the legacy rescalePseudoReaction shim, which
%   dispatches the lipid special case.
%
% Usage: model = scaleBioMass(model, component, new_value, balance_out, dispOutput)

if nargin < 5
    dispOutput = true;
end
if nargin < 4
    balance_out = '';
end

% Current fractions (uses the shared yeast biomass config).
[X, P, C, R, D, L, I, F] = sumBioMass(model, false);
content_all = {'biomass','carbohydrate','protein','lipid','RNA','DNA','ion','cofactor'};
content_Cap = {'X','C','P','L','R','D','I','F'};
pos       = strcmp(content_all, component);
old_value = eval(content_Cap{pos});
f         = new_value / old_value;
model     = rescalePseudoReaction(model, component, f);

if ~isempty(balance_out)
    X             = sumBioMass(model, false);
    pos           = strcmp(content_all, balance_out);
    balance_value = eval(content_Cap{pos});
    f             = (balance_value + (1 - X)) / balance_value;
    model         = rescalePseudoReaction(model, balance_out, f);
end
sumBioMass(model, dispOutput);
end

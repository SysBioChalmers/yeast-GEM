function model = addSBOterms(model)
% addSBOterms  yeast-GEM shim — delegates to RAVEN's assignSBOterms.
%
%   The legacy implementation had a typo: the pseudoreaction-SBO loop
%   used `for i = numel(model.rxns)` (single iteration over the last
%   reaction) instead of `1:numel(model.rxns)`. To keep this function
%   byte-equivalent to the pre-refactor output (and avoid spurious
%   SBO churn in saveYeastModel diffs), we pass
%   `onlyLastReactionForPseudo = true`. Flip the flag off here if
%   yeast-GEM ever decides to start tagging every pseudoreaction.
%
% Usage: model = addSBOterms(model)

model = assignSBOterms(model, struct('onlyLastReactionForPseudo', true));
end

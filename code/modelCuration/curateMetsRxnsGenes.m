function newModel = curateMetsRxnsGenes(model, metsInfo, genesInfo, rxnsCoeffs, rxnsInfo)
% curateMetsRxnsGenes  Yeast-GEM batch curation entry point.
%
%   Thin wrapper around RAVEN's ``curateModelFromTables`` that pins the
%   yeast-GEM id prefixes (``'s_'`` for new metabolites, ``'r_'`` for
%   new reactions). All four TSV file arguments and the matching
%   semantics are identical to the generic upstream function — see
%   ``curateModelFromTables`` for the full docstring.
%
%   Existing v8_*/v9_* curation scripts (and TEMPLATEcuration) call
%   this function with 1–4 arguments; the shim preserves the
%   original signature so no caller needs to change.
%
%   Requires RAVEN ≥ the commit that added core/curateModelFromTables.m
%   (currently the feat/yeast-gem-shared branch).
%
%   Input:
%       model       RAVEN model structure to be curated.
%       metsInfo    relative path to the *.tsv file with metabolite
%                   information, or 'none' to skip.
%       genesInfo   relative path to the *.tsv file with gene
%                   information, or 'none'.
%       rxnsCoeffs  relative path to the *.tsv file with reaction
%                   stoichiometric coefficients, or 'none'.
%       rxnsInfo    relative path to the *.tsv file with reaction
%                   information, or 'none'.
%
%   Output:
%       newModel    curated RAVEN model structure.
%
% Usage: newModel = curateMetsRxnsGenes(model, metsInfo, genesInfo, ...
%                       rxnsCoeffs, rxnsInfo)

if nargin == 4
    error('Provide both a ''rxnsInfo'' and a ''rxnsCoeffs'' file')
end
if nargin < 5
    rxnsInfo = 'none';
end
if nargin < 4
    rxnsCoeffs = 'none';
end
if nargin < 3
    genesInfo = 'none';
end

newModel = curateModelFromTables(model, metsInfo, genesInfo, ...
    rxnsCoeffs, rxnsInfo, 's_', 'r_');
end

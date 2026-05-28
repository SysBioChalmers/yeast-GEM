function regenerate()
% regenerate
%   Produce the reference bundle for the Python-vs-MATLAB CI gates.
%   Run from this directory with RAVEN on the path:
%
%       cd code/python/tests/reference
%       regenerate
%
%   Writes (in this directory):
%       yeast-GEM.xml       MATLAB-produced canonical model
%       metrics.yml         reference values for the level-2 gate
%       provenance.yml      MATLAB / RAVEN / solver / git SHA
%
%   This script is intentionally a thin orchestrator. The behaviours it
%   captures (load, commit pipeline, growth/essentialGenes metrics) are
%   already implemented under code/. See PORTING_PLAN.md, validation
%   strategy, for what the bundle is used for.

% TODO: implement once the Python port reaches the level-2 gate. For
% phase 1 of the port (scaffold), this file documents the contract:
%
%   1. Load model with loadYeastModel.
%   2. Run the canonical commit pipeline (saveYeastModel /
%      commitYeastModel) into yeast-GEM.xml in this folder.
%   3. Run growth, essentialGenes, anaerobic_flux_predictions,
%      sumBioMass and persist the numeric results to metrics.yml.
%   4. Record MATLAB / RAVEN / solver / git SHA to provenance.yml.
error('regenerate.m is a phase-1 scaffold; not implemented yet.');
end

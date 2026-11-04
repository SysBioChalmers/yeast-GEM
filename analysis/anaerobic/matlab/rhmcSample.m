function [X, diag] = rhmcSample(model, nSamples, opts)
% rhmcSample  Uniform flux samples with constrained Riemannian HMC.
%   [X, diag] = rhmcSample(model, nSamples, opts)
%
%   Samples {v : S*v = 0, lb <= v <= ub} uniformly with PolytopeSamplerMatlab
%   (Kook et al. 2022, NeurIPS; the sampler behind COBRA's 'RHMC' option).
%   opts.cap limits internal (non-exchange) reaction bounds to [-cap, cap]
%   (the +-1000 defaults only give loop directions an arbitrary length);
%   opts.seed, opts.maxTime are passed to the sampler.
%
%   X     nRxns-by-nDraws samples
%   diag  table per reaction: split-R-hat and effective sample size over 4
%         consecutive segments of the chain, plus run time.

if nargin < 3, opts = struct(); end
if ~isfield(opts, 'cap'),     opts.cap = inf;    end
if ~isfield(opts, 'seed'),    opts.seed = 1;     end
if ~isfield(opts, 'maxTime'), opts.maxTime = 36000; end

[~, exch] = getExchangeRxns(model);
internal = true(numel(model.rxns), 1); internal(exch) = false;
lb = model.lb; ub = model.ub;
lb(internal) = max(lb(internal), -opts.cap);
ub(internal) = min(ub(internal),  opts.cap);

P = struct('Aeq', sparse(model.S), 'beq', zeros(size(model.S,1),1), 'lb', lb, 'ub', ub);
so = default_options();
so.seed = opts.seed;
so.maxTime = opts.maxTime;
tic;
o = sample(P, nSamples, so);
t = toc;
X = o.samples;

nd = 4*floor(size(X,2)/4);
rhat = zeros(size(X,1),1); ess = rhat;
for i = 1:size(X,1)
    [rhat(i), ess(i)] = mcmcDiag(reshape(X(i,1:nd), [], 4));
end
diag = table(model.rxns, rhat, ess, min(X,[],2), max(X,[],2), 'VariableNames', {'rxn','rhat','ess','min','max'});
diag.Properties.UserData = struct('seconds', t, 'draws', size(X,2));
end

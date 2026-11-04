function [samples, info] = chrrChains(model, nSamples, thinning, nBurnin, nChains, seed, opts)
% chrrChains  Coordinate hit-and-run with rounding, several chains from one setup.
%   [samples, info] = chrrChains(model, nSamples, thinning, nBurnin, nChains, seed)
%
%   The algorithm of RAVEN v3.0.0-rc1 sampleCHRR (Haraldsdottir et al. 2017,
%   Bioinformatics 33:1741): fold stoichiometrically fixed reactions into the
%   equalities, reduce to the null space, round with the maximum-volume
%   inscribed ellipsoid, then walk with coordinate hit-and-run. The setup
%   (flux variability and rounding) is done once and shared by all chains,
%   which start from the same interior point with different seeds.
%
%   Two additions to RAVEN's version, for strongly anisotropic polytopes:
%   opts.cap limits the bounds of internal (non-exchange) reactions to
%   [-cap, cap] (the arbitrary +-1000 bounds make loop directions 1e5 times
%   longer than, e.g., the growth range), and the rounding is repeated on the
%   rounded polytope, as in COBRA's CHRR preprocessing, until the last
%   ellipsoid is close to the unit ball (opts.maxRounds, opts.maxIter).
%   opts.fixedTol (default 1e-7, RAVEN's value) is the flux-range width below
%   which a reaction is treated as fixed.
%
%   samples  nRxns-by-nSamples-by-nChains
%   info     .nDimensions, .mveConverged, .roundCond, .fixedRxns, .tSetup, .tSteps

if nargin < 7, opts = struct(); end
if ~isfield(opts, 'cap'),       opts.cap = inf;   end
if ~isfield(opts, 'maxRounds'), opts.maxRounds = 8; end
if ~isfield(opts, 'maxIter'),   opts.maxIter = 400; end
if ~isfield(opts, 'fixedTol'),  opts.fixedTol = 1e-7; end
tic;
nRxns = numel(model.rxns);
[~, exch] = getExchangeRxns(model);
internal = true(nRxns, 1); internal(exch) = false;
model.lb(internal) = max(model.lb(internal), -opts.cap);
model.ub(internal) = min(model.ub(internal),  opts.cap);
S  = full(model.S);
lb = model.lb(:);
ub = model.ub(:);

[minF, maxF] = getAllowedBounds(model);
minF = minF(:); maxF = maxF(:);
fixed = (maxF - minF) < opts.fixedTol;
fixedIdx = find(fixed);
Aeq = [S; sparse(1:numel(fixedIdx), fixedIdx, 1, numel(fixedIdx), nRxns)];
beq = [zeros(size(S,1),1); (maxF(fixedIdx) + minF(fixedIdx)) / 2];
v0 = lsqminnorm(full(Aeq), beq);
N  = null(full(Aeq));
d  = size(N, 2);

free = ~fixed;
Nf = N(free, :);
A_full = [Nf; -Nf];
b_full = [ub(free) - v0(free); v0(free) - lb(free)];
keep = sqrt(sum(A_full.^2, 2)) > 1e-10;
A_full = A_full(keep, :);
b_full = b_full(keep);

x0 = sampleChebyshevCenter(A_full, b_full);
% Repeated rounding: round, transform, and round the result again until the
% new ellipsoid is nearly the unit ball (condition number of its matrix < 3).
center = zeros(d, 1); E = eye(d);
A_r = A_full; b_r = b_full; start = x0;
roundCond = [];
for r = 1:opts.maxRounds
    [c_k, E_k, converged] = sampleMaxVolEllipse(A_r, b_r, start, opts.maxIter);
    center = center + E*c_k;
    E = E*E_k;
    A_r = A_full * E;
    b_r = b_full - A_full*center;
    % A non-converged ellipsoid's centre need not be strictly interior; walking
    % from a point outside the polytope makes almost every hit-and-run step
    % empty. Re-centre on the Chebyshev centre of the rounded polytope.
    y0 = sampleChebyshevCenter(A_r, b_r);
    center = center + E*y0;
    b_r = b_full - A_full*center;
    roundCond(end+1) = cond(E_k); %#ok<AGROW>
    fprintf('  rounding %d: converged %d, cond(E_k) %.3g, min slack at centre %.3g\n', ...
        r, converged, roundCond(end), min(b_r));
    start = zeros(d, 1);
    if converged && roundCond(end) < 3, break; end
end

assert(all(b_r > 0), 'chrrChains: walk would start outside the polytope');
info.nDimensions = d;
info.mveConverged = converged;
info.roundCond = roundCond;
info.fixedRxns = model.rxns(fixedIdx);
info.tSetup = toc;

tic;
samples = zeros(nRxns, nSamples, nChains);
skipped = zeros(nChains, 1);
parfor c = 1:nChains
    rs = RandStream('mt19937ar', 'Seed', seed + c);
    [samples(:,:,c), skipped(c)] = walk(A_r, b_r, center, E, N, v0, d, nSamples, thinning, nBurnin, rs);
end
info.tSteps = toc;
info.skippedFraction = skipped / (nBurnin + thinning*nSamples);
end

function [sol, skipped] = walk(A_r, b_r, center, E, N, v0, d, nSamples, thinning, nBurnin, rs)
y = zeros(d, 1);
s = b_r - A_r*y;
sol = zeros(numel(v0), nSamples);
rec = 0;
skipped = 0;
colTol = 1e-12;
for stepno = 1:(nBurnin + thinning*nSamples)
    j = randi(rs, d);
    col = A_r(:, j);
    pos = col > colTol;
    neg = col < -colTol;
    ratios = s ./ col;
    if any(pos); tHi = min(ratios(pos)); else; tHi = inf;  end
    if any(neg); tLo = max(ratios(neg)); else; tLo = -inf; end
    if ~isfinite(tHi) || ~isfinite(tLo) || tHi <= tLo
        skipped = skipped + 1;
        continue;
    end
    t = tLo + rand(rs)*(tHi - tLo);
    y(j) = y(j) + t;
    s = s - col*t;
    if stepno > nBurnin && mod(stepno - nBurnin, thinning) == 0
        rec = rec + 1;
        sol(:, rec) = v0 + N*(center + E*y);
    end
end
end

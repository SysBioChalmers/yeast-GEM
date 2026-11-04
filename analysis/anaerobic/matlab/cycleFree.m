function W = cycleFree(model, V, keepIdx)
% cycleFree  Remove internal cycles from sampled flux distributions.
%   W = cycleFree(model, V, keepIdx)
%
%   For each column v of V (reactions-by-samples), solves
%       min  sum_i |w_i|
%       s.t. S*w = 0,
%            w_k = v_k               for k in keepIdx (exchanges, biomass),
%            0 <= w_i <= v_i         where v_i >= 0,
%            v_i <= w_i <= 0         where v_i < 0,
%   the CycleFreeFlux projection (Desouki et al. 2015, Bioinformatics
%   31:2159). W carries the same exchange fluxes as V with every internal
%   cycle removed, and never reverses or enlarges a flux. Uses Gurobi
%   directly; columns are processed in parallel if a pool is open.

S = model.S;
n = size(S, 2);
keep = false(n, 1); keep(keepIdx) = true;
W = zeros(size(V));
base.A = sparse(S);
base.rhs = zeros(size(S, 1), 1);
base.sense = '=';
base.modelsense = 'min';
params.OutputFlag = 0;
params.Threads = 1;
params.FeasibilityTol = 1e-9;
tol = 1e-7;
parfor k = 1:size(V, 2)
    v = V(:, k);
    g = base;
    % Sampled fluxes balance only to ~1e-10, so the bounds get a small
    % tolerance; exact bounds could make the projection spuriously infeasible.
    g.lb = min(v, 0) - tol;
    g.ub = max(v, 0) + tol;
    g.lb(keep) = v(keep) - tol;
    g.ub(keep) = v(keep) + tol;
    g.obj = sign(v) .* ~keep;           % |w_i| = sign(v_i)*w_i given the sign bounds
    r = gurobi(g, params);
    if ~strcmp(r.status, 'OPTIMAL')
        error('cycleFree: projection of sample %d returned %s', k, r.status);
    end
    W(:, k) = r.x;
end
end

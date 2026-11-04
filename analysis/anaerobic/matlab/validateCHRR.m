function validateCHRR()
% validateCHRR  Check RAVEN's sampleCHRR against polytopes with a known
%   uniform distribution.
%
%   Simplex tests: reactions R1..Rn draw a metabolite that a fixed reaction
%   supplies at rate 1, with weights w_i, so sum(w_i*v_i) = 1 and v_i >= 0.
%   The uniform distribution on that simplex gives w_i*v_i ~ Beta(1, n-1).
%   Unequal weights make the polytope elongated, which is what the rounding
%   step has to cope with. Each case runs 4 chains with different seeds, and
%   reports the error in mean and variance, the Kolmogorov-Smirnov distance
%   to Beta(1, n-1), and the potential scale reduction R-hat.

paths = setupStudy('3.0.0-rc1');
out = fullfile(paths.results, 'chrr_validation');
if ~isfolder(out), mkdir(out); end

cases = {
    'simplex3',          [1 1 1]
    'simplex3_elongated',[1 1e-3 1]
    'simplex20',         ones(1,20)
    'simplex20_mixed',   [ones(1,10) 1e-2*ones(1,5) 1e2*ones(1,5)]
    };
nSamples = 2000; thinning = 200; nChains = 4;
res = table();
for c = 1:size(cases,1)
    w = cases{c,2}; n = numel(w);
    model = simplexModel(w);
    chains = zeros(n, nSamples, nChains);
    for k = 1:nChains
        [sol, info] = sampleCHRR(model, nSamples, thinning, 1000, k);
        chains(:,:,k) = sol(2:end,:) .* w(:);       % scaled so each ~ Beta(1,n-1)
    end
    pooled = reshape(chains, n, []);
    a = 1; b = n-1;
    muTrue = a/(a+b); varTrue = a*b/((a+b)^2*(a+b+1));
    for i = 1:n
        x = pooled(i,:);
        [~,~,ks] = kstest(x, 'CDF', makedist('Beta','a',a,'b',b));
        res = [res; table(cases(c,1), i, n, mean(x), muTrue, var(x), varTrue, ks, ...
            rhat(squeeze(chains(i,:,:))), info.mveConverged, ...
            'VariableNames', {'case','rxn','n','mean','mean_true','var','var_true','ks','rhat','mve_converged'})]; %#ok<AGROW>
    end
    fprintf('%s: max |mean err| %.4f (true %.4f), max KS %.4f, max Rhat %.4f\n', cases{c,1}, ...
        max(abs(res.mean(strcmp(res.case,cases{c,1})) - muTrue)), muTrue, ...
        max(res.ks(strcmp(res.case,cases{c,1}))), max(res.rhat(strcmp(res.case,cases{c,1}))));
end
writetable(res, fullfile(out, 'chrr_validation.tsv'), 'FileType','text','Delimiter','\t');
end

function model = simplexModel(w)
n = numel(w);
model.id = 'simplex'; model.name = 'simplex';
model.rxns = [{'R0'}; strcat('R', cellstr(num2str((1:n)','%d')))];
model.rxnNames = model.rxns;
model.mets = {'A'}; model.metNames = {'A'};
model.comps = {'c'}; model.compNames = {'c'}; model.metComps = 1;
model.S = sparse([1, -w]);
model.lb = [1; zeros(n,1)]; model.ub = [1; 1000*ones(n,1)];
model.rev = zeros(n+1,1); model.c = zeros(n+1,1); model.b = 0;
model.genes = {}; model.grRules = repmat({''}, n+1, 1); model.rxnGeneMat = sparse(n+1,0);
end

function R = rhat(x)
% Gelman-Rubin potential scale reduction for one quantity, x is draws-by-chains.
[n, m] = size(x);
B = n*var(mean(x,1)); W = mean(var(x,0,1));
R = sqrt(((n-1)/n*W + B/n)/W);
end

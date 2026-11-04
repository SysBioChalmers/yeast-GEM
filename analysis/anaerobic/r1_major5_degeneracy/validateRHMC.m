function validateRHMC()
% validateRHMC  Check PolytopeSamplerMatlab (constrained Riemannian HMC,
%   Kook et al. 2022) against polytopes with a known uniform distribution:
%   the simplex cases sum(w_i*v_i) = 1, v_i >= 0, where
%   the uniform distribution gives w_i*v_i ~ Beta(1, n-1).

paths = setupStudy('3.0.0-rc1');
addRHMC();
out = fullfile(paths.results, 'rhmc_validation');
if ~isfolder(out), mkdir(out); end

cases = {
    'simplex3',          [1 1 1]
    'simplex3_elongated',[1 1e-3 1]
    'simplex20',         ones(1,20)
    'simplex20_mixed',   [ones(1,10) 1e-2*ones(1,5) 1e2*ones(1,5)]
    'simplex200_mixed',  [ones(1,100) 1e-3*ones(1,50) 1e3*ones(1,50)]
    };
res = table();
for c = 1:size(cases,1)
    w = cases{c,2}; n = numel(w);
    P = struct('Aeq', w, 'beq', 1, 'lb', zeros(n,1), 'ub', inf(n,1));
    opts = default_options(); opts.seed = c; opts.maxTime = 3600;
    tic; o = sample(P, 4000, opts); t = toc;
    x = o.samples .* w(:);                       % each ~ Beta(1, n-1)
    a = 1; b = n-1; muTrue = a/(a+b); varTrue = a*b/((a+b)^2*(a+b+1));
    nd = 4*floor(size(x,2)/4);          % one chain, split into 4 segments for R-hat/ESS
    for i = 1:n
        [~,~,ks] = kstest(x(i,:), 'CDF', makedist('Beta','a',a,'b',b));
        [rh, es] = mcmcDiag(reshape(x(i,1:nd), [], 4));
        res = [res; table(cases(c,1), i, n, mean(x(i,:)), muTrue, var(x(i,:)), varTrue, ks, ...
            es, rh, size(x,2), t, ...
            'VariableNames', {'case','var','n','mean','mean_true','var_sample','var_true','ks','ess','rhat','draws','seconds'})]; %#ok<AGROW>
    end
    sel = strcmp(res.case, cases{c,1});
    fprintf('%-18s n=%3d  draws %6d  %5.0f s | mean rel err max %.3f  var rel err max %.3f  KS max %.4f  ESS min %.0f  Rhat max %.3f\n', ...
        cases{c,1}, n, size(x,2), t, max(abs(res.mean(sel)-muTrue))/muTrue, max(abs(res.var_sample(sel)-varTrue))/varTrue, ...
        max(res.ks(sel)), min(res.ess(sel)), max(res.rhat(sel)));
end
writetable(res, fullfile(out, 'rhmc_validation.tsv'), 'FileType','text','Delimiter','\t');
end

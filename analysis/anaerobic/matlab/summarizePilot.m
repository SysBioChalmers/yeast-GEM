function summarizePilot(file)
% summarizePilot  Distributions of key fluxes in a pilot sample, against the
%   pFBA solution of the same model.

paths = setupStudy('3.0.0-rc1');
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
model = anaerobicModel(loadTaggedModel(paths.yeast, '9.1.0'));
model = setParam(model, 'eq', 'r_1714', -23);
sol = solveLP(model, 1);
S = load(file, 'W');
keyRxns = {'r_4041','growth'; 'r_1808','glycerol exch'; 'r_1761','ethanol exch'; 'r_1672','CO2 exch'; ...
           'r_2056','succinate exch'; 'r_1634','acetate exch'; 'r_0713','Mdh1 (mito)'; 'r_0714','Mdh2 (cyto, blocked)'; ...
           'r_0658','Idp1 (mito)'; 'r_0659','Idp2 (cyto, blocked)'; 'r_0471','Gdh1/3 (NADP)'; 'r_0470','Gdh2 (NAD)'; ...
           'r_0476','Gln1'; 'r_0472','Glt1 (blocked)'; 'r_0226','ATP synthase (mito)'; 'r_0227','PM H+-ATPase'; ...
           'r_1654','NH4 exch'; 'r_0491','Gpd1/2'};
fprintf('%-24s %10s | %10s %10s %10s %10s | %s\n', 'reaction', 'pFBA', 'q2.5', 'median', 'q97.5', 'P(|v|>1e-6)', 'pFBA inside 95%?');
for k = 1:size(keyRxns,1)
    i = find(strcmp(model.rxns, keyRxns{k,1}));
    v = S.W(i,:); q = prctile(v, [2.5 50 97.5]);
    fprintf('%-24s %10.4f | %10.4f %10.4f %10.4f %10.2f | %d\n', keyRxns{k,2}, sol.x(i), q, mean(abs(v) > 1e-6), sol.x(i) >= q(1) && sol.x(i) <= q(3));
end
end

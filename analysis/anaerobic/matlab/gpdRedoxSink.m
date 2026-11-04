function gpdRedoxSink()
% gpdRedoxSink  How v9.1.0 + new constraints reoxidises NADH without
%   glycerol formation: compare pFBA wild type with gpd1 gpd2 (glucose 19.1),
%   list the largest changes in reactions that use NADH or NAD+, and all
%   exchange fluxes that change.

paths = setupStudy('3.0.0-rc1');
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
m = anaerobicModel(loadTaggedModel(paths.yeast, '9.1.0'));
m = setParam(m, 'eq', 'r_1714', -19.1);
m = setParam(m, 'obj', 'r_4041', 1);
wt = solveLP(m, 1);
ko = solveLP(setParam(removeGenes(m, {'YDL022W','YOL059W'}), 'eq', 'r_1714', -19.1), 1);
nadh = find(ismember(m.metNames, {'NADH','NAD'}));
usesNad = any(m.S(nadh,:) ~= 0, 1)';
d = ko.x - wt.x;
[~, o] = sort(abs(d) .* usesNad, 'descend');
fprintf('Largest changes in NAD(H) reactions, gpd1 gpd2 minus wild type:\n');
for k = o(1:20)'
    eq = constructEquations(m, m.rxns(k));
    fprintf('%-8s %-45s wt %8.3f ko %8.3f | %s\n', m.rxns{k}, m.rxnNames{k}(1:min(45,end)), wt.x(k), ko.x(k), eq{1});
end
[~, exch] = getExchangeRxns(m);
fprintf('\nExchange fluxes that change by more than 0.01:\n');
for k = exch(:)'
    if abs(d(k)) > 0.01
        fprintf('%-8s %-35s wt %8.3f ko %8.3f\n', m.rxns{k}, m.rxnNames{k}(1:min(35,end)), wt.x(k), ko.x(k));
    end
end
end

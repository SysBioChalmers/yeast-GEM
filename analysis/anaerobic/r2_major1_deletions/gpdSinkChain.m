function gpdSinkChain()
% gpdSinkChain  In v9.1.0 + new constraints with gpd1 gpd2 deleted (glucose
%   19.1), repeatedly close the new product exchange that pFBA uses in place
%   of glycerol, and report growth and the next one, until growth stops or
%   eight exits are closed. Maps the chain of alternative NADH sinks.

paths = setupStudy('3.0.0-rc1');
log = startLog(paths, 'gpdSinkChain'); %#ok<NASGU>
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
m = applyEnvironment(loadTaggedModel(paths.yeast, 'dev'), 'anaerobic');
m = setParam(m, 'eq', 'r_1714', -19.1);
m = setParam(m, 'obj', 'r_4041', 1);
wt = solveLP(m, 1);
[~, exch] = getExchangeRxns(m);
% Products of the wild type (secreted, > 1e-6) stay allowed.
wtOut = exch(wt.x(exch) > 1e-6);
ko = removeGenes(m, {'YDL022W','YOL059W'});
closed = {};
for it = 1:8
    sol = solveLP(ko, 1);
    if isempty(sol.x), fprintf('step %d: infeasible\n', it); break; end
    mu = sol.x(strcmp(ko.rxns, 'r_4041'));
    newOut = setdiff(exch(sol.x(exch) > 1e-3), wtOut);
    [~, j] = max(sol.x(newOut));
    if mu < 1e-6 || isempty(newOut)
        fprintf('step %d: growth %.4f, no new product exchange\n', it, mu); break
    end
    r = newOut(j);
    fprintf('step %d: growth %.4f; largest new product %s (%s) %.3f\n', it, mu, ko.rxns{r}, ko.rxnNames{r}, sol.x(r));
    closed{end+1} = ko.rxns{r}; %#ok<AGROW>
    ko.ub(r) = 0;
end
fprintf('closed: %s\n', strjoin(closed, ' '));
sol = solveLP(ko, 1);
if ~isempty(sol.x)
    fprintf('exchanges that differ from the wild type by > 0.01 after closing:\n');
    for k = exch(:)'
        if abs(sol.x(k) - wt.x(k)) > 0.01
            fprintf('  %-8s %-30s wt %8.3f ko %8.3f\n', ko.rxns{k}, ko.rxnNames{k}(1:min(30,end)), wt.x(k), sol.x(k));
        end
    end
end
end

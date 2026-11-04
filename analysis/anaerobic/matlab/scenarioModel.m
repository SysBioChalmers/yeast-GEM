function [model, spec] = scenarioModel(paths, variant, scenario, optFrac)
% scenarioModel  Build the constrained model for one (variant, scenario).
%   variant   see variantModel, e.g. 'v910_new', 'v902_old', 'v911c_new'
%   scenario  'sjoberg'             anaerobic batch, glucose 23 (Figure 1)
%             'bakker_wt'           anaerobic batch, glucose 19.1 (Bakker 2000, Table 5)
%             'bakker_adh3'         same, adh3 deletion, glucose 12.6
%             'bakker_adh3adh4'     same, adh3 adh4 deletion (in silico control)
%             'boles_wt'            anaerobic glucose-limited chemostat, D = 0.10
%             'boles_mae1'          same, mae1 deletion
%             'flux_<dataset>'      a dataset of flux_data_anaerobic.tsv, glucose as measured
%
%   Batch and fluxomics scenarios fix glucose uptake and constrain growth to
%   its optimum; chemostat scenarios fix growth and constrain glucose uptake to
%   its minimum. optFrac (default 0.99999) sets how close: at >= 0.9999 the
%   objective is fixed (the optimal face, the set of optimal solutions);
%   below that it is a bound (e.g. 0.99 for 'at least 99% of the maximum').
%   optFrac > 1 sets a band of +-(optFrac-1) around both reference values
%   instead: e.g. 1.01 keeps glucose uptake within 99-101% of the measured
%   rate and growth within 99-101% of its maximum at the measured rate (batch),
%   or growth within 99-101% of D and glucose within 99-101% of its minimum
%   (chemostat). spec records what was set.

if nargin < 4, optFrac = 0.99999; end
model = variantModel(paths, variant);

spec = struct('variant', variant, 'scenario', scenario, 'ko', {{}}, 'optFrac', optFrac, ...
    'mode', 'batch', 'glucose', NaN, 'growth', NaN);
switch scenario
    case 'sjoberg',         spec.glucose = 23;
    case 'bakker_wt',       spec.glucose = 19.1;
    case 'bakker_adh3',     spec.glucose = 12.6; spec.ko = {'YMR083W'};
    case 'bakker_adh3adh4', spec.glucose = 12.6; spec.ko = {'YMR083W','YGL256W'};
    case 'boles_wt',        spec.mode = 'chemostat'; spec.growth = 0.10;
    case 'boles_mae1',      spec.mode = 'chemostat'; spec.growth = 0.10; spec.ko = {'YKL029C'};
    otherwise
        assert(startsWith(scenario, 'flux_'), 'Unknown scenario %s', scenario);
        T = readtable(fullfile(paths.yeast, 'data', 'physiology', 'flux_data_anaerobic.tsv'), ...
            'FileType', 'text', 'ReadVariableNames', false, 'Delimiter', '\t');
        row = strcmp(T{:,6}, extractAfter(scenario, 'flux_')) & strcmp(T{:,7}, 'r_1714');
        assert(nnz(row) == 1, 'No single glucose row for %s', scenario);
        spec.glucose = T{row, 4};
end

if ~isempty(spec.ko)
    missing = setdiff(spec.ko, model.genes);
    assert(isempty(missing), 'Genes not in model: %s', strjoin(missing, ' '));
    model = removeGenes(model, spec.ko);
end

bio = 'r_4041';
if strcmp(spec.mode, 'batch')
    model = setParam(model, 'eq', 'r_1714', -spec.glucose);
    model = setParam(model, 'obj', bio, 1);
    sol = solveLP(model, 1);
    if isempty(sol.x)                       % e.g. a lethal knockout
        spec.growthMax = 0; spec.pfba = []; spec.feasible = false;
        return
    end
    spec.growthMax = sol.x(getIndexes(model, bio, 'rxns'));
    spec.feasible = spec.growthMax > 1e-6;  % no growth: nothing near-optimal to sample
    spec.pfba = sol.x;
    if optFrac > 1
        d = optFrac - 1;
        model = setParam(model, 'lb', 'r_1714', -spec.glucose * (1 + d));
        model = setParam(model, 'ub', 'r_1714', -spec.glucose * (1 - d));
        model = setParam(model, 'lb', bio, (1 - d) * spec.growthMax);
        model = setParam(model, 'ub', bio, (1 + d) * spec.growthMax);
    elseif optFrac >= 0.9999
        model = setParam(model, 'eq', bio, optFrac * spec.growthMax);
    else
        model = setParam(model, 'lb', bio, optFrac * spec.growthMax);
    end
else
    % Chemostat: growth fixed at D; the least glucose that supports it.
    model = setParam(model, 'lb', 'r_1714', -1000);
    model = setParam(model, 'eq', bio, spec.growth);
    model = setParam(model, 'obj', 'r_1714', 1);     % maximise the (negative) exchange
    sol = solveLP(model, 1);
    if isempty(sol.x)
        spec.glucoseMin = NaN; spec.pfba = []; spec.feasible = false;
        return
    end
    spec.feasible = true;
    spec.glucoseMin = -sol.x(getIndexes(model, 'r_1714', 'rxns'));
    spec.pfba = sol.x;
    if optFrac > 1
        d = optFrac - 1;
        model = setParam(model, 'lb', bio, spec.growth * (1 - d));
        model = setParam(model, 'ub', bio, spec.growth * (1 + d));
        model = setParam(model, 'lb', 'r_1714', -spec.glucoseMin * (1 + d));
        model = setParam(model, 'ub', 'r_1714', -spec.glucoseMin * (1 - d));
    elseif optFrac >= 0.9999
        model = setParam(model, 'eq', 'r_1714', -spec.glucoseMin / optFrac);
    else
        model = setParam(model, 'lb', 'r_1714', -spec.glucoseMin / optFrac);
    end
end
model.c(:) = 0;
end

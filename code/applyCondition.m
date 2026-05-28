function model = applyCondition(model, name)
% applyCondition  Apply a named condition preset to the yeast-GEM model.
%
%   model = applyCondition(model, name) loads data/conditions/<name>.yml
%   and applies the described modifications (bound changes, biomass
%   stoichiometry deltas, amino acid ratio switches, ...) to model.
%
%   Available presets:
%       'minimal_Y6'         minimal media (replaces minimal_Y6.m)
%       'anaerobic'          anaerobic conditions (replaces anaerobicModel.m)
%       'glycine_nitrogen'   glycine as sole N source (replaces glycineNitrogenSource.m)
%       'nitrogen_limitation' N-limited (replaces nitrogenLimitation.m)
%
%   This is the data-as-code refactor described in PORTING_PLAN.md
%   phase 2. The YAML files in data/conditions/ are the single source
%   of truth shared with the Python loader
%   (yeastgem.conditions.apply).
%
% Usage: model = applyCondition(model, 'anaerobic')

funcDir = fileparts(mfilename('fullpath'));
yamlPath = fullfile(funcDir, '..', 'data', 'conditions', [name '.yml']);
if ~isfile(yamlPath)
    error('applyCondition:unknownCondition', ...
        'No such condition: %s (looked for %s)', name, yamlPath);
end
cond = readYAML(yamlPath);

% --- Step 1: prelude ----------------------------------------------------
if isfield(cond, 'prelude') && isfield(cond.prelude, 'reset_exchanges')
    [~, exchangeRxns] = getExchangeRxns(model, cond.prelude.reset_exchanges);
    model.lb(exchangeRxns) = 0;
    model.ub(exchangeRxns) = 1000;
end

% --- Step 2: cofactor pseudoreaction edits ------------------------------
if isfield(cond, 'cofactor_pseudoreaction')
    cp = cond.cofactor_pseudoreaction;
    cofacIdx = getIndexes(model, cp.rxn_id, 'rxns');
    if isfield(cp, 'remove_mets')
        for i = 1:numel(cp.remove_mets)
            metId = cp.remove_mets{i}.met;
            metIdx = getIndexes(model, metId, 'mets');
            model.S(metIdx, cofacIdx) = 0;
        end
    end
    if isfield(cp, 'charge_balance_met')
        balanceIdx = find(strcmp(model.mets, cp.charge_balance_met));
        model.S(balanceIdx, cofacIdx) = 0;
        model.S(balanceIdx, cofacIdx) = ...
            -sum(model.S(:, cofacIdx) .* model.metCharges, 'omitnan');
    end
end

% --- Step 3: amino acid ratio switch ------------------------------------
if isfield(cond, 'amino_acid_ratio')
    aerobic = strcmp(cond.amino_acid_ratio, 'aerobic');
    model = changeAminoAcidRatio(model, aerobic);
end

% --- Step 4: biomass stoichiometry delta --------------------------------
if isfield(cond, 'biomass_stoichiometry_delta')
    delta = cond.biomass_stoichiometry_delta;
    bioIdx = getIndexes(model, delta.rxn_id, 'rxns');
    if isfield(delta, 'add')
        for i = 1:numel(delta.add)
            entry = delta.add{i};
            metIdx = getIndexes(model, entry.met, 'mets');
            model.S(metIdx, bioIdx) = full(model.S(metIdx, bioIdx)) + entry.coef;
        end
    end
end

% --- Step 5: bounds -----------------------------------------------------
nUptake = 0;
if isfield(cond, 'bounds')
    for i = 1:numel(cond.bounds)
        b = cond.bounds{i};
        rxnIdx = find(strcmp(model.rxns, b.rxn));
        if isempty(rxnIdx)
            warning('applyCondition:missingRxn', ...
                'Reaction %s not found in model; skipping.', b.rxn);
            continue;
        end
        if isfield(b, 'lb')
            model.lb(rxnIdx) = b.lb;
            if b.lb < 0
                nUptake = nUptake + 1;
            end
        end
        if isfield(b, 'ub')
            model.ub(rxnIdx) = b.ub;
        end
    end
end

% --- Step 6: sanity check on uptake count (mirrors minimal_Y6 warning) --
if isfield(cond, 'expected_uptake_count')
    if nUptake ~= cond.expected_uptake_count
        warning('applyCondition:uptakeCountMismatch', ...
            'Expected %d uptake reactions, applied %d. Some may be missing from the model.', ...
            cond.expected_uptake_count, nUptake);
    end
end

end

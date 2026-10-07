function model = applyEnvironment(model,environment)
% applyEnvironment
%   Applies an environment (medium, oxygen availability, nitrogen source,
%   ...) to yeast-GEM, as defined in data/conditions/<environment>.yml.
%   The same files are read by the Python yeastgem.conditions.apply, so
%   both languages apply identical constraints. The function needs no
%   other toolbox: the files are read with a small reader that supports
%   the subset of YAML used in data/conditions.
%
%   Steps, in this order (as yeastgem.conditions.apply):
%   amino_acid_ratio             'aerobic' or 'anaerobic' column of
%                                data/physiology/aminoAcid_Bjorkeroth2020.tsv
%   prelude.reset_exchanges      'out', 'in' or 'all': set those exchange
%                                reactions to lb 0, ub 1000
%   cofactor_pseudoreaction      remove metabolites from the pseudoreaction
%                                and recompute the charge-balancing H+
%   biomass_stoichiometry_delta  add coefficients to a reaction
%   bounds                       set lb and/or ub of listed reactions
%   expected_uptake_count        warn if fewer lb -1000 bounds were applied
%
% Input:
%   model           yeast-GEM model structure
%   environment     name of a file in data/conditions, without .yml
%                   (e.g. 'anaerobic', 'minimal_Y6', 'glycine_nitrogen',
%                   'nitrogen_limitation'), or the path to
%                   such a file
%
% Output:
%   model           model with the environment applied
%
% Usage: model = applyEnvironment(model,environment)

codeDir = fileparts(mfilename('fullpath'));
if isfile(environment)
    envFile = environment;
else
    envFile = fullfile(codeDir,'..','data','conditions',[environment '.yml']);
    if ~isfile(envFile)
        error('applyEnvironment:unknownEnvironment', ...
            'No environment file %s.', envFile);
    end
end
env = readEnvironmentFile(envFile);

% Amino acid ratio of the protein pseudoreaction
if isfield(env,'amino_acid_ratio')
    switch env.amino_acid_ratio
        case 'aerobic',   aerobic = true;
        case 'anaerobic', aerobic = false;
        otherwise
            error('applyEnvironment:aminoAcidRatio', ...
                'amino_acid_ratio must be aerobic or anaerobic.');
    end
    oldPath = addpath(fullfile(codeDir,'otherChanges'));
    restorePath = onCleanup(@() path(oldPath));
    model = changeAminoAcidRatio(model,aerobic);
    clear restorePath
end

% Reset exchange reactions (exchange rule of RAVEN's getExchangeRxns)
if isfield(env,'prelude') && isfield(env.prelude,'reset_exchanges')
    noProducts   = ~any(model.S > 0,1)';
    noSubstrates = ~any(model.S < 0,1)';
    switch lower(env.prelude.reset_exchanges)
        case 'out', exch = noProducts;
        case 'in',  exch = noSubstrates;
        otherwise,  exch = noProducts | noSubstrates;
    end
    model.lb(exch) = 0;
    model.ub(exch) = 1000;
end

% Cofactor pseudoreaction: remove metabolites, recompute the H+ balance
if isfield(env,'cofactor_pseudoreaction')
    cp = env.cofactor_pseudoreaction;
    rxnIdx = findIndex(model.rxns,cp.rxn_id);
    if isfield(cp,'remove_mets')
        for i = 1:numel(cp.remove_mets)
            model.S(findIndex(model.mets,cp.remove_mets{i}.met),rxnIdx) = 0;
        end
    end
    if isfield(cp,'charge_balance_met')
        balIdx = findIndex(model.mets,cp.charge_balance_met);
        model.S(balIdx,rxnIdx) = 0;
        metIdx = find(model.S(:,rxnIdx));
        unknown = metIdx(isnan(model.metCharges(metIdx)));
        if ~isempty(unknown)
            error('applyEnvironment:unknownCharge', ...
                'Cannot charge balance %s, no charge for: %s', ...
                cp.rxn_id, strjoin(model.mets(unknown),', '));
        end
        model.S(balIdx,rxnIdx) = ...
            -sum(full(model.S(metIdx,rxnIdx)).*model.metCharges(metIdx));
    end
end

% Add coefficients to the biomass (or another) reaction
if isfield(env,'biomass_stoichiometry_delta')
    delta = env.biomass_stoichiometry_delta;
    rxnIdx = findIndex(model.rxns,delta.rxn_id);
    for i = 1:numel(delta.add)
        metIdx = findIndex(model.mets,delta.add{i}.met);
        model.S(metIdx,rxnIdx) = model.S(metIdx,rxnIdx) + delta.add{i}.coef;
    end
end

% Reaction bounds
nUptake = 0;
if isfield(env,'bounds')
    for i = 1:numel(env.bounds)
        b = env.bounds{i};
        rxnIdx = find(strcmp(model.rxns,b.rxn));
        if isempty(rxnIdx)
            warning('applyEnvironment:missingRxn', ...
                'Reaction %s not in the model; skipped.', b.rxn);
            continue
        end
        if isfield(b,'lb')
            model.lb(rxnIdx) = b.lb;
            nUptake = nUptake + (b.lb == -1000);
        end
        if isfield(b,'ub')
            model.ub(rxnIdx) = b.ub;
        end
    end
end
if isfield(env,'expected_uptake_count') && nUptake ~= env.expected_uptake_count
    warning('applyEnvironment:uptakeCount', ...
        'Expected %d uptake reactions, applied %d.', env.expected_uptake_count, nUptake);
end
end

function idx = findIndex(list,id)
idx = find(strcmp(list,id));
if isempty(idx)
    error('applyEnvironment:missingId','%s not in the model.',id);
end
end

%% Reader for the YAML subset used in data/conditions
% Supported: block mappings, block sequences of scalars or flow mappings
% ({ key: value, ... }), quoted and plain scalars, numbers, true/false/null
% and # comments. Anything else is an error rather than a silent misread.

function data = readEnvironmentFile(fileName)
rawLines = regexp(fileread(fileName),'\r?\n','split');
lines = {}; indents = [];
for i = 1:numel(rawLines)
    content = stripComment(rawLines{i});
    if isempty(strtrim(content)) || any(strcmp(strtrim(content),{'---','...'}))
        continue
    end
    if contains(content,sprintf('\t'))
        error('applyEnvironment:yaml','%s: tabs are not allowed (line %d).',fileName,i);
    end
    lines{end+1} = strtrim(content); %#ok<AGROW>
    indents(end+1) = numel(content) - numel(strtrim([content 'x'])) + 1; %#ok<AGROW>
end
if isempty(lines)
    data = struct();
    return
end
[data,k] = parseBlock(lines,indents,1,indents(1));
if k <= numel(lines)
    error('applyEnvironment:yaml','Unexpected indentation at "%s".',lines{k});
end
end

function [value,k] = parseBlock(lines,indents,k,ind)
if startsWith(lines{k},'-')
    value = {};
    while k <= numel(lines) && indents(k) == ind && startsWith(lines{k},'-')
        item = strtrim(lines{k}(2:end));
        if isempty(item)
            error('applyEnvironment:yaml','Nested block in sequence not supported.');
        end
        value{end+1,1} = parseInline(item); %#ok<AGROW>
        k = k + 1;
    end
else
    value = struct();
    while k <= numel(lines) && indents(k) == ind && ~startsWith(lines{k},'-')
        tok = splitKey(lines{k});
        if isempty(tok)
            error('applyEnvironment:yaml','Cannot read "%s".',lines{k});
        end
        key = tok{1};
        if isempty(tok{2})
            if k < numel(lines) && (indents(k+1) > ind || ...
                    (indents(k+1) == ind && startsWith(lines{k+1},'-')))
                [value.(key),k] = parseBlock(lines,indents,k+1,indents(k+1));
            else
                value.(key) = [];
                k = k + 1;
            end
        else
            value.(key) = parseInline(tok{2});
            k = k + 1;
        end
    end
end
if k <= numel(lines) && indents(k) > ind
    error('applyEnvironment:yaml','Unexpected indentation at "%s".',lines{k});
end
end

function value = parseInline(s)
s = strtrim(s);
if startsWith(s,'{') && endsWith(s,'}')
    value = struct();
    parts = splitTopLevel(s(2:end-1));
    for i = 1:numel(parts)
        tok = splitKey(parts{i});
        if isempty(tok)
            error('applyEnvironment:yaml', ...
                'Cannot read "%s" in "%s" (quote values that contain commas).',parts{i},s);
        end
        value.(tok{1}) = parseScalar(tok{2});
    end
elseif startsWith(s,'[') && endsWith(s,']')
    parts = splitTopLevel(s(2:end-1));
    value = cellfun(@parseScalar,parts(:),'UniformOutput',false);
else
    value = parseScalar(s);
end
end

function tok = splitKey(s)
% {key, value} for 'key: value' or 'key:', empty if s is not a key.
tok = regexp(s,'^([A-Za-z_]\w*):(.*)$','tokens','once');
if ~isempty(tok)
    if isempty(tok{2})
        tok{2} = '';
    elseif ~isspace(tok{2}(1))   % 'key:value' is not a key in YAML
        tok = {};
    else
        tok{2} = strtrim(tok{2});
    end
end
end

function value = parseScalar(s)
s = strtrim(s);
if numel(s) >= 2 && s(1) == '"' && s(end) == '"'
    value = regexprep(s(2:end-1),'\\(["\\])','$1');
elseif numel(s) >= 2 && s(1) == '''' && s(end) == ''''
    value = strrep(s(2:end-1),'''''','''');
elseif ~isempty(regexp(s,'^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$','once'))
    value = str2double(s);
elseif any(strcmp(s,{'true','True','TRUE'}))
    value = true;
elseif any(strcmp(s,{'false','False','FALSE'}))
    value = false;
elseif isempty(s) || any(strcmp(s,{'null','Null','NULL','~'}))
    value = [];
else
    value = s;
end
end

function parts = splitTopLevel(s)
parts = {}; depth = 0; quote = ''; start = 1;
for i = 1:numel(s)
    c = s(i);
    if ~isempty(quote)
        if c == quote, quote = ''; end
    elseif c == '"' || c == ''''
        quote = c;
    elseif c == '{' || c == '['
        depth = depth + 1;
    elseif c == '}' || c == ']'
        depth = depth - 1;
    elseif c == ',' && depth == 0
        parts{end+1} = strtrim(s(start:i-1)); %#ok<AGROW>
        start = i + 1;
    end
end
last = strtrim(s(start:end));
if ~isempty(last)
    parts{end+1} = last;
end
end

function line = stripComment(line)
quote = '';
for i = 1:numel(line)
    c = line(i);
    if ~isempty(quote)
        if c == quote, quote = ''; end
    elseif c == '"' || c == ''''
        quote = c;
    elseif c == '#' && (i == 1 || any(line(i-1) == [' ' sprintf('\t')]))
        line = line(1:i-1);
        return
    end
end
end

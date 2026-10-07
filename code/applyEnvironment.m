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
%                                (a metabolite without charge counts as 0,
%                                with a warning)
%   biomass_stoichiometry_delta  add coefficients to a reaction
%   bounds                       set lb and/or ub of listed reactions
%   expected_uptake_count        warn if fewer lb -1000 bounds were applied
%
%   The environment is checked before the model is changed: an invalid
%   value (e.g. a non-numeric bound, lb > ub, an unknown reset_exchanges
%   value) is an error and nothing is applied.
%
% Input:
%   model           yeast-GEM model structure
%   environment     name of a file in data/conditions, without .yml
%                   (e.g. 'anaerobic', 'minimal_Y6', 'glycine_nitrogen',
%                   'nitrogen_limitation', 'carnitine'), or the path to
%                   such a file
%
% Output:
%   model           model with the environment applied
%
% Usage: model = applyEnvironment(model,environment)

environment = char(environment);
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
if isfield(model,'unconstrained') && any(model.unconstrained)
    error('applyEnvironment:boundaryMets', ...
        'Models with boundary metabolites (model.unconstrained) are not supported.');
end

% Check the environment before changing the model
aerobic = [];
if isfield(env,'amino_acid_ratio') && ~isempty(env.amino_acid_ratio)
    if ~ischar(env.amino_acid_ratio) || ~any(strcmp(env.amino_acid_ratio,{'aerobic','anaerobic'}))
        error('applyEnvironment:aminoAcidRatio', ...
            'amino_acid_ratio must be aerobic or anaerobic.');
    end
    aerobic = strcmp(env.amino_acid_ratio,'aerobic');
end
exch = false(numel(model.rxns),1);
if isfield(env,'prelude') && isfield(env.prelude,'reset_exchanges') ...
        && ~isempty(env.prelude.reset_exchanges)
    reset = env.prelude.reset_exchanges;
    if ~ischar(reset) || ~any(strcmpi(reset,{'in','out','all'}))
        error('applyEnvironment:resetExchanges', ...
            'reset_exchanges must be in, out or all.');
    end
    % Exchange rule of RAVEN's getExchangeRxns
    noProducts   = ~any(model.S > 0,1)';
    noSubstrates = ~any(model.S < 0,1)';
    switch lower(reset)
        case 'out', exch = noProducts;
        case 'in',  exch = noSubstrates;
        otherwise,  exch = noProducts | noSubstrates;
    end
end
newBounds = checkBounds(model,env,exch);

% Amino acid ratio of the protein pseudoreaction
if ~isempty(aerobic)
    oldPath = addpath(fullfile(codeDir,'otherChanges'));
    restorePath = onCleanup(@() path(oldPath));
    model = changeAminoAcidRatio(model,aerobic);
    clear restorePath
end

% Reset exchange reactions
model.lb(exch) = 0;
model.ub(exch) = 1000;

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
            warning('applyEnvironment:unknownCharge', ...
                'Charge balance of %s: no charge for %s, counted as 0.', ...
                cp.rxn_id, strjoin(model.mets(unknown),', '));
        end
        model.S(balIdx,rxnIdx) = ...
            -sum(full(model.S(metIdx,rxnIdx)).*model.metCharges(metIdx),'omitnan');
    end
end

% Add coefficients to the biomass (or another) reaction
if isfield(env,'biomass_stoichiometry_delta')
    delta = env.biomass_stoichiometry_delta;
    rxnIdx = findIndex(model.rxns,delta.rxn_id);
    if isfield(delta,'add')
        for i = 1:numel(delta.add)
            metIdx = findIndex(model.mets,delta.add{i}.met);
            model.S(metIdx,rxnIdx) = model.S(metIdx,rxnIdx) + delta.add{i}.coef;
        end
    end
end

% Reaction bounds
model.lb(newBounds.idx) = newBounds.lb;
model.ub(newBounds.idx) = newBounds.ub;
nUptake = sum(newBounds.hasLb & newBounds.lb == -1000);
if isfield(env,'expected_uptake_count') && nUptake ~= env.expected_uptake_count
    warning('applyEnvironment:uptakeCount', ...
        'Expected %d uptake reactions, applied %d.', env.expected_uptake_count, nUptake);
end
end

function nb = checkBounds(model,env,exch)
% The bounds to set (reaction index, lb, ub, whether lb was given), checked:
% bounds must be finite numbers and lb <= ub, given the exchange reset.
nb = struct('idx',zeros(0,1),'lb',zeros(0,1),'ub',zeros(0,1),'hasLb',false(0,1));
if ~isfield(env,'bounds')
    return
end
for i = 1:numel(env.bounds)
    b = env.bounds{i};
    rxnIdx = find(strcmp(model.rxns,b.rxn));
    if isempty(rxnIdx)
        warning('applyEnvironment:missingRxn', ...
            'Reaction %s not in the model; skipped.', b.rxn);
        continue
    end
    if exch(rxnIdx), lb = 0; ub = 1000; else, lb = model.lb(rxnIdx); ub = model.ub(rxnIdx); end
    for key = {'lb','ub'}
        if isfield(b,key{1})
            v = b.(key{1});
            if ~isnumeric(v) || ~isscalar(v) || ~isreal(v) || ~isfinite(v)
                error('applyEnvironment:bound', ...
                    '%s: %s must be a finite number.', b.rxn, key{1});
            end
        end
    end
    if isfield(b,'lb'), lb = b.lb; end
    if isfield(b,'ub'), ub = b.ub; end
    if lb > ub
        error('applyEnvironment:bound','%s: lb %g is larger than ub %g.', b.rxn, lb, ub);
    end
    nb.idx(end+1,1) = rxnIdx;
    nb.lb(end+1,1) = lb;
    nb.ub(end+1,1) = ub;
    nb.hasLb(end+1,1) = isfield(b,'lb');
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
tok = regexp(s,'^([A-Za-z]\w*):(.*)$','tokens','once');
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
parts = {}; depth = 0; quote = ''; start = 1; i = 1;
while i <= numel(s)
    c = s(i);
    if ~isempty(quote)
        i = quoteStep(s,i,quote);
        if i < 0, quote = ''; i = -i; end
    else
        if (c == '"' || c == '''') && quoteStarts(s,i)
            quote = c;
        elseif c == '{' || c == '['
            depth = depth + 1;
        elseif c == '}' || c == ']'
            depth = depth - 1;
        elseif c == ',' && depth == 0
            parts{end+1} = strtrim(s(start:i-1)); %#ok<AGROW>
            start = i + 1;
        end
        i = i + 1;
    end
end
last = strtrim(s(start:end));
if ~isempty(last)
    parts{end+1} = last;
end
end

function line = stripComment(line)
quote = ''; i = 1;
while i <= numel(line)
    c = line(i);
    if ~isempty(quote)
        i = quoteStep(line,i,quote);
        if i < 0, quote = ''; i = -i; end
    else
        if (c == '"' || c == '''') && quoteStarts(line,i)
            quote = c;
        elseif c == '#' && (i == 1 || any(line(i-1) == [' ' sprintf('\t')]))
            line = line(1:i-1);
            return
        end
        i = i + 1;
    end
end
end

function i = quoteStep(s,i,quote)
% Next position inside a quoted scalar; negative if the quote closes at i.
% Escapes: \x in double quotes, '' in single quotes.
if quote == '"' && s(i) == '\'
    i = i + 2;
elseif s(i) == quote && quote == '''' && i < numel(s) && s(i+1) == ''''
    i = i + 2;
elseif s(i) == quote
    i = -(i + 1);
else
    i = i + 1;
end
end

function tf = quoteStarts(s,i)
% A quote character starts a quoted scalar only at the start of a value
% (after '{', '[', ',', ': ' or a sequence '- '); elsewhere, as in
% 5'-phosphate, it is part of a plain scalar.
j = i - 1;
while j >= 1 && isspace(s(j))
    j = j - 1;
end
tf = j < 1 || any(s(j) == '{[,:') || (s(j) == '-' && (j == 1 || isspace(s(j-1))));
end

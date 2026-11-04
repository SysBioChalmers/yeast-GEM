function testApplyEnvironment
% testApplyEnvironment  Error cases, the combination of bound entries and the
%   YAML reader of applyEnvironment, on a small model. Errors if a check fails.
%   Run with the code folder on the path:
%   matlab -batch "addpath('code'); addpath('code/python/tests/reference'); testApplyEnvironment"

model = toyModel();
tmp = tempname; mkdir(tmp);
cleanup = onCleanup(@() rmdir(tmp,'s'));
envFile = @(text) writeFile([tempname(tmp) '.yml'],text);

% reset_exchanges selects by the side without products or substrates
expected = struct('out',[1 0 0 0],'in',[0 1 0 0],'all',[1 1 0 0]);
for name = fieldnames(expected)'
    m = applyEnvironment(model,envFile(sprintf('prelude:\n  reset_exchanges: %s\n',name{1})));
    changed = (m.lb ~= model.lb | m.ub ~= model.ub)';
    assert(isequal(changed,logical(expected.(name{1}))),'reset_exchanges %s',name{1});
    assert(all(m.lb(changed) == 0) && all(m.ub(changed) == 1000));
end

% an identifier that is not in the model is an error
expectError(@() applyEnvironment(model,envFile(sprintf('bounds:\n  - { rxn: nope, lb: 0 }\n'))), ...
    'applyEnvironment:missingId');
expectError(@() applyEnvironment(model,envFile(sprintf( ...
    'cofactor_pseudoreaction:\n  rxn_id: P\n  remove_mets:\n    - { met: nope }\n'))), ...
    'applyEnvironment:missingId');
expectError(@() applyEnvironment(model,envFile(sprintf( ...
    'biomass_stoichiometry_delta:\n  rxn_id: P\n  add:\n    - { met: nope, coef: 1 }\n'))), ...
    'applyEnvironment:missingId');

% entries for one reaction are combined in order
m = applyEnvironment(model,envFile(sprintf( ...
    'bounds:\n  - { rxn: R, lb: -1 }\n  - { rxn: R, ub: 3 }\n  - { rxn: R, lb: -2 }\n')));
assert(m.lb(3) == -2 && m.ub(3) == 3,'combined bounds');

% lb above ub, a non-numeric bound
expectError(@() applyEnvironment(model,envFile(sprintf('bounds:\n  - { rxn: R, lb: 5, ub: 1 }\n'))), ...
    'applyEnvironment:bound');
expectError(@() applyEnvironment(model,envFile(sprintf('bounds:\n  - { rxn: R, ub: high }\n'))), ...
    'applyEnvironment:bound');

% the uptake count warns only on a difference
cfg = sprintf('expected_uptake_count: %d\nbounds:\n  - { rxn: R, lb: -1000 }\n',1);
lastwarn(''); applyEnvironment(model,envFile(cfg)); [~,id] = lastwarn;
assert(isempty(id),'no warning when the count is as expected');
lastwarn(''); applyEnvironment(model,envFile(strrep(cfg,'count: 1','count: 2'))); [~,id] = lastwarn;
assert(strcmp(id,'applyEnvironment:uptakeCount'));

% a metabolite without charge counts as 0, with a warning
lastwarn('');
m = applyEnvironment(model,envFile(sprintf( ...
    'cofactor_pseudoreaction:\n  rxn_id: P\n  charge_balance_met: h\n')));
[~,id] = lastwarn;
assert(strcmp(id,'applyEnvironment:unknownCharge'));
assert(abs(full(m.S(4,4)) - 1) < 1e-12,'charge balance');

% numbers are read as decimal numbers: 1e3 is 1000, 010 is 10
m = applyEnvironment(model,envFile(sprintf('bounds:\n  - { rxn: R, lb: 010, ub: 1e3 }\n')));
assert(m.lb(3) == 10 && m.ub(3) == 1000,'numbers');

% a mapping in a sequence must be a flow mapping
expectError(@() applyEnvironment(model,envFile(sprintf('bounds:\n  - rxn: R\n    lb: 0\n'))), ...
    'applyEnvironment:yaml');
disp('testApplyEnvironment: all checks passed');
end

function model = toyModel()
% A -> B with an exchange of A (A =>), an import of B (=> B) and H+.
model.mets = {'a';'b';'x';'h'};
model.rxns = {'EX_a';'IN_b';'R';'P'};
model.S = sparse([-1 0 -1 -1; 0 1 1 1; 0 0 0 -1; 0 0 0 0.5]);
model.lb = [-5;-3;-10;0];
model.ub = [7;4;10;1000];
model.metCharges = [0;-1;NaN;1];
end

function file = writeFile(file,text)
fid = fopen(file,'w'); fprintf(fid,'%s',text); fclose(fid);
end

function expectError(fn,id)
try
    fn();
catch err
    assert(strcmp(err.identifier,id),'expected %s, got %s',id,err.identifier);
    return
end
error('expected the error %s',id);
end

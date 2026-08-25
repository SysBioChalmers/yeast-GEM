function comparePhase2(preDir, postDir)
% comparePhase2  Compare pre- vs post-refactor model structs for phase 2.
%
%   Loads <preDir>/<name>.mat and <postDir>/<name>.mat for each of the
%   four conditions and checks equality of the fields the phase-2
%   refactor touches (lb, ub, S, metCharges). Prints a summary; exits
%   the MATLAB session with code 0 if all match, 1 otherwise.

warning('off','all');
restoredefaultpath; rehash toolboxcache;

conds = {'minimal_Y6', 'anaerobicModel', 'glycineNitrogenSource', 'nitrogenLimitation'};
tol = 1e-12;
allOk = true;

for i = 1:numel(conds)
    name = conds{i};
    fprintf('=== %s ===\n', name);
    pre = load(fullfile(preDir, [name '.mat']));
    post = load(fullfile(postDir, [name '.mat']));
    P = pre.model; Q = post.model;

    okRxns = isequal(P.rxns, Q.rxns);
    okMets = isequal(P.mets, Q.mets);
    okLb   = isequal(P.lb, Q.lb);
    okUb   = isequal(P.ub, Q.ub);
    okS    = nnz(abs(P.S - Q.S) > tol) == 0;

    fprintf('  rxns:  %s\n', tf(okRxns));
    fprintf('  mets:  %s\n', tf(okMets));
    fprintf('  lb:    %s', tf(okLb));
    if ~okLb
        diffs = find(P.lb ~= Q.lb);
        fprintf(' (%d rxn(s) differ: %s)', numel(diffs), ...
            strjoin(P.rxns(diffs(1:min(end,5)))', ', '));
    end
    fprintf('\n');
    fprintf('  ub:    %s', tf(okUb));
    if ~okUb
        diffs = find(P.ub ~= Q.ub);
        fprintf(' (%d rxn(s) differ: %s)', numel(diffs), ...
            strjoin(P.rxns(diffs(1:min(end,5)))', ', '));
    end
    fprintf('\n');
    fprintf('  S:     %s', tf(okS));
    if ~okS
        [iMet, iRxn] = find(abs(P.S - Q.S) > tol);
        fprintf(' (%d entry/entries differ; first: met=%s rxn=%s pre=%g post=%g)', ...
            numel(iMet), P.mets{iMet(1)}, P.rxns{iRxn(1)}, ...
            full(P.S(iMet(1), iRxn(1))), full(Q.S(iMet(1), iRxn(1))));
    end
    fprintf('\n');

    if ~(okRxns && okMets && okLb && okUb && okS)
        allOk = false;
    end
end

fprintf('\n');
if allOk
    fprintf('OVERALL: all four conditions semantically equal (pre vs post).\n');
    exit(0);
else
    fprintf('OVERALL: MISMATCH detected. See above.\n');
    exit(1);
end
end

function s = tf(b)
if b
    s = 'OK';
else
    s = 'DIFFER';
end
end

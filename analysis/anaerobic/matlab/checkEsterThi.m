function checkEsterThi()
% checkEsterThi  Two evidence-based direction checks on v9.1.0:
%   E  Iah1 (YOR126C) ethyl-ester esterases hydrolysis only (no reverse
%      ester synthesis): effect on the growth benchmark, in particular the
%      N-limited anaerobic chemostats that export ethyl pyruvate.
%   T  Thi5-family HMP-P synthase (r_1603) needs O2: blocked in the
%      anaerobic model; does anaerobic growth then depend on thiamine uptake?

paths = setupStudy('3.0.0-rc1');
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
base = loadTaggedModel(paths.yeast, '9.1.0');
est = {'r_4711','r_4712','r_4715','r_4716','r_4717','r_4720','r_4721'};
eqs = constructEquations(base, est);
for k = 1:numel(est), fprintf('%s [%g,%g] %s\n', est{k}, base.lb(strcmp(base.rxns,est{k})), base.ub(strcmp(base.rxns,est{k})), eqs{k}); end
E = setParam(base, 'ub', est, 0);                  % written as ester synthesis; keep hydrolysis only

grp = [repmat({'Nlim_aer'},9,1); repmat({'Clim_aer'},11,1); repmat({'Clim_anaer'},6,1); repmat({'Nlim_anaer'},6,1)];
for mdl = {'released', base; 'E', E}'
    figure('Visible','off'); [~, e, g] = growthNewSim(mdl{2}); close all
    o = double(e(:,4)); pr = g(:,4);
    fprintf('\n%s: growth R2 %.4f MRE %.4f\n', mdl{1}, 1 - sum((o-pr).^2)/sum((o-mean(o)).^2), mean(abs(pr-o)./o));
    for gname = unique(grp)'
        sel = strcmp(grp, gname{1});
        fprintf('  %-11s MRE %.4f\n', gname{1}, mean(abs(pr(sel)-o(sel))./o(sel)));
    end
    sel = strcmp(grp, 'Nlim_anaer');
    fprintf('  N-lim anaerobic measured %s\n  predicted %s\n', sprintf('%.3f ', o(sel)), sprintf('%.3f ', pr(sel)));
end

an = setParam(setParam(anaerobicModel(base), 'eq', 'r_1714', -23), 'obj', 'r_4041', 1);
[~, exch] = getExchangeRxns(an);
thi = exch(contains(lower(an.rxnNames(exch)), {'thiamin','hmp','pyrimidine'}));
for k = thi', fprintf('\n%s %s [%g,%g]', an.rxns{k}, an.rxnNames{k}, an.lb(k), an.ub(k)); end
fprintf('\n');
s = solveLP(an); fprintf('anaerobic Sjoberg growth %.4f; r_1603 flux %.2e\n', s.x(strcmp(an.rxns,'r_4041')), s.x(strcmp(an.rxns,'r_1603')));
s = solveLP(setParam(an, 'ub', 'r_1603', 0)); g = 0; if ~isempty(s.x), g = s.x(strcmp(an.rxns,'r_4041')); end
fprintf('with r_1603 blocked: growth %.4f\n', g);
end

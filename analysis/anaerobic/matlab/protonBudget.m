function protonBudget(variant)
% protonBudget  Cytosolic H+ budget of the anaerobic Sjoberg simulation (pFBA,
%   glucose 23): net H+[c] production per reaction, grouped as glycolysis to
%   ethanol, glycerol, ammonium uptake and assimilation, biomass, transport over
%   the plasma membrane, mitochondrial transport and other. Prints the groups
%   and the largest single contributions; also the budget of the lumped
%   overall reactions glucose -> 2 ethanol + 2 CO2 and glucose -> glycerol.
if nargin < 1, variant = 'v910_new'; end
paths = setupStudy('3.0.0-rc1');
cd(fullfile(paths.yeast,'code','otherChanges'));
m = variantModel(paths, variant);
m = setParam(m,'eq','r_1714',-23);
s = solveLP(m,1);
h = strcmp(m.mets,'s_0794');                       % H+[c]
contrib = full(m.S(h,:))' .* s.x;                  % mmol H+[c]/gDW/h produced (+) or consumed (-)
ss = cellfun(@joinSS, m.subSystems, 'UniformOutput', false);
fprintf('%s, growth %.4f\n', variant, s.x(strcmp(m.rxns,'r_4041')));
grp = repmat({'other'}, numel(m.rxns), 1);
grp(contains(ss,'Glycolysis') | contains(ss,'Pyruvate metabolism')) = {'glycolysis/pyruvate'};
grp(contains(ss,'Glycerolipid') | ismember(m.rxns,{'r_0491','r_0489'})) = {'glycerol'};
grp(ismember(m.rxns,{'r_0227'})) = {'plasma membrane H+-ATPase'};
grp(ismember(m.rxns,{'r_1115','r_0470','r_0471','r_0476'})) = {'NH4 uptake and assimilation'};
grp(startsWith(m.rxnNames,'growth') | ismember(m.rxns,{'r_4041','r_4047','r_4048','r_4049','r_4050','r_4598','r_4599'})) = {'biomass pseudoreactions'};
isPM = any(m.S(strcmp(m.mets,'s_0796'),:),1)';   % involves H+[e]
grp(isPM & ~strcmp(grp,'plasma membrane H+-ATPase')) = {'other H+ transport over plasma membrane'};
isMito = any(m.S(strcmp(m.mets,'s_0799'),:),1)'; % involves H+[m]
grp(isMito & strcmp(grp,'other')) = {'exchange with mitochondria'};
[g,~,gi] = unique(grp);
tot = accumarray(gi, contrib);
for i = 1:numel(g), fprintf('  %-42s %+8.3f\n', g{i}, tot(i)); end
fprintf('  %-42s %+8.3f\n', 'sum (should be 0)', sum(contrib));
[~,o] = sort(abs(contrib),'descend');
fprintf('Largest contributions:\n');
for i = o(1:15)'
    fprintf('  %-8s %-55s flux %8.3f  H+[c] %+8.3f\n', m.rxns{i}, m.rxnNames{i}, s.x(i), contrib(i));
end
end

function s = joinSS(x)
if iscell(x), s = strjoin(x(:)','; '); else, s = char(x); end
end

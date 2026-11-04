function c1Mechanism()
% c1Mechanism  Network inventory for the bound-free curation options
%   (v9.1.0 + new anaerobic constraints, Sjoberg batch, glucose 23):
%   1. every reaction of the one-carbon, glycine, threonine, propanol,
%      glycolaldehyde and ethyl pyruvate genes/metabolites, with bounds,
%      equation and GPR;
%   2. pFBA flux of those reactions with glycine cleavage (GCV) blocked (as
%      released), GCV open in both directions as encoded, and GCV forward only;
%   3. the flux range of each GCV reaction at optimal growth with GCV open.

paths = setupStudy('3.0.0-rc1');
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
m = anaerobicModel(loadTaggedModel(paths.yeast, '9.1.0'));
m = setParam(setParam(m, 'eq', 'r_1714', -23), 'obj', 'r_4041', 1);
gcv = {'r_0501','r_0507','r_0509'};

genes = {'YDR019C','YMR189W','YAL044C','YFL018C', ...          % GCV1 GCV2 GCV3 LPD1
         'YBR263W','YLR058C','YEL046C','YFL030W', ...           % SHM1 SHM2 GLY1 AGX1
         'YGR204W','YBR084W','YKR080W','YGL125W','YPL023C', ... % ADE3 MIS1 MTD1 MET13 MET12
         'YOR388C','YPL262W','YJR139C','YER086W','YCL064C', ...  % FDH1 FUM1 HOM6 ILV1 CHA1
         'YOR375C'};                                              % GDH1
metNames = {'L-glycine','glycolaldehyde','formate','propanal','1-propanol','ethyl pyruvate', ...
            '5,10-methylenetetrahydrofolate','10-formyl-THF','2-oxobutanoate','L-threonine'};
idx = find(ismember(m.rxns, gcv) | any(m.rxnGeneMat(:, ismember(m.genes, genes)), 2) | ...
           any(m.S(ismember(m.metNames, metNames), :) ~= 0, 1)');

variants = {'released (GCV blocked)', m};
mo = setParam(m, 'ub', gcv, 1000);
mo.lb(ismember(mo.rxns, gcv)) = -1000 * (m.rev(ismember(m.rxns, gcv)) == 1);
variants(end+1, :) = {'GCV open as encoded', mo};
mf = setParam(setParam(m, 'ub', gcv, 1000), 'lb', gcv, 0);
variants(end+1, :) = {'GCV forward only', mf};

X = zeros(numel(idx), size(variants,1));
for v = 1:size(variants,1)
    s = solveLP(variants{v,2}, 1);
    X(:, v) = s.x(idx);
    fprintf('%-24s growth %.4f  glycerol %.3f  glycine ex %.4f  glycolaldehyde ex %.4f  formate ex %.4f  CO2 %.2f\n', ...
        variants{v,1}, s.x(strcmp(m.rxns,'r_4041')), s.x(strcmp(m.rxns,'r_1808')), s.x(strcmp(m.rxns,'r_1810')), ...
        s.x(strcmp(m.rxns,'r_1814')), s.x(strcmp(m.rxns,'r_1665')), s.x(strcmp(m.rxns,'r_1672')));
end
fprintf('\nrxn      lb     ub    | flux released / GCV as encoded / GCV forward | name | GPR | equation\n');
eqs = constructEquations(m, m.rxns(idx));
for k = 1:numel(idx)
    i = idx(k);
    fprintf('%-8s %6g %6g | %8.4f %8.4f %8.4f | %s | %s | %s\n', m.rxns{i}, m.lb(i), m.ub(i), X(k,:), ...
        m.rxnNames{i}, m.grRules{i}, eqs{k});
end

fprintf('\nGCV flux range at optimal growth, GCV open as encoded:\n');
s = solveLP(mo); mu = s.x(strcmp(mo.rxns,'r_4041'));
mm = setParam(mo, 'lb', 'r_4041', 0.99999*mu);
for r = gcv
    lo = solveLP(setParam(mm, 'obj', r{1}, -1)); hi = solveLP(setParam(mm, 'obj', r{1}, 1));
    fprintf('  %s [%.4f, %.4f]\n', r{1}, lo.x(strcmp(mm.rxns,r{1})), hi.x(strcmp(mm.rxns,r{1})));
end
end

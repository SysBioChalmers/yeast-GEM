function kozakTest()
% kozakTest  Kozak et al. 2014: acs1 acs2 + E. coli PFL (IMZ383), anaerobic
%   glucose-limited chemostat, D 0.10. Measured biomass yield 0.082 vs 0.096
%   g/g for the reference; the authors attribute the difference to weak-acid
%   uncoupling by the secreted formate (0.18) and acetate (0.11
%   mmol/gDW/h). Here: biomass yield (growth fixed at 0.10, glucose uptake
%   minimised) of the reference and of the PFL strain, the latter with free
%   and with forced formate and acetate secretion.
paths = setupStudy('3.0.0-rc1');
cd(fullfile(paths.yeast,'code','otherChanges'));
an = variantModel(paths,'v910_new');
pfl = removeGenes(an, {'YAL054C','YLR153C'}, false, false, false);
r = struct('rxns',{{'PFL'}}, 'rxnNames',{{'pyruvate formate lyase'}}, ...
    'equations',{{'pyruvate[c] + coenzyme A[c] => acetyl-CoA[c] + formate[c]'}});
pfl = addRxns(pfl, r, 3);
forced = setParam(pfl, 'lb', {'r_1793','r_1634'}, [0.18 0.11]);
for c = {{'reference', an}, {'acs1 acs2 + PFL, free', pfl}, {'acs1 acs2 + PFL, formate 0.18 and acetate 0.11 forced', forced}}
    m = setParam(c{1}{2}, 'eq', 'r_4041', 0.10);
    m = setParam(m, 'lb', 'r_1714', -1000);
    m = setParam(m, 'obj', 'r_1714', 1);
    s = solveLP(m, 1);
    if isempty(s.x), fprintf('%-58s infeasible\n', c{1}{1}); continue; end
    q = -s.x(strcmp(m.rxns,'r_1714'));
    fprintf('%-58s glucose %.3f, biomass yield %.4f g/g, formate %.3f, acetate %.3f, ATPase %.3f\n', c{1}{1}, q, 0.10/(q*0.18016), ...
        s.x(strcmp(m.rxns,'r_1793')), s.x(strcmp(m.rxns,'r_1634')), s.x(strcmp(m.rxns,'r_0227')));
end
fprintf('measured: reference 0.096 +- 0.002, acs1 acs2 + PFL 0.082 +- 0.002 g/g\n');
end

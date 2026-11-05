function frdKnockouts()
% frdKnockouts  Anaerobic growth (Sjoberg batch) of frd1, osm1 and frd1 osm1 in
%   v9.1.1 with the new anaerobic constraints; frd1 osm1 cannot grow
%   anaerobically, single mutants can (Arikawa et al. 1998).
paths = setupStudy('3.0.0-rc1');
log = startLog(paths, 'frdKnockouts'); %#ok<NASGU>
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
an = variantModel(paths, 'v911_new');
for ko = {{}, {'YEL047C'}, {'YJR051W'}, {'YEL047C','YJR051W'}}
    m = an; if ~isempty(ko{1}), m = removeGenes(m, ko{1}); end
    m = setParam(setParam(m, 'eq', 'r_1714', -23), 'obj', 'r_4041', 1);
    s = solveLP(m, 1); g = 0; if ~isempty(s.x), g = s.x(strcmp(m.rxns,'r_4041')); end
    fprintf('knockout %-16s growth %.4f\n', strjoin(ko{1}, '+'), g);
end
end

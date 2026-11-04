function nlimExports()
% nlimExports  Exports used by v9.1.0 in the N-limited anaerobic chemostats of
%   the growth benchmark (same set-up as growth.m: anaerobicModel, protein
%   0.28, RNA 0.0329 balanced by carbohydrate, glutamate synthase and glycine
%   cleavage reopened), listing those outside the known-product list.

paths = setupStudy('3.0.0-rc1');
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
m = anaerobicModel(loadTaggedModel(paths.yeast, '9.1.0'));
m = scaleBioMass(m, 'protein', 0.28, '', false);
m = scaleBioMass(m, 'RNA', 0.0329, 'carbohydrate', false);
m = setParam(m, 'ub', {'r_0472','r_0501','r_0507','r_0509'}, 1000);
products = {'r_1761','r_1672','r_1808','r_1634','r_2056','r_2033','r_1546','r_1551','r_2100','r_1832','r_4041','r_2111'};
fid = fopen(fullfile(paths.yeast,'data','physiology','chemostatData_Tobias2013.tsv'));
d = textscan(fid, '%f32 %f32 %f32 %f32', 'Delimiter', '\t', 'HeaderLines', 1); fclose(fid);
d = double([d{:}]); d = d(27:32, :);            % N-limited anaerobic rows, as in growth.m
[~, exch] = getExchangeRxns(m);
for i = 1:size(d,1)
    mt = m;
    pos = {'r_1714','r_1992','r_1654'};
    for j = 1:3
        if abs(d(i,j)) == 1000, mt = setParam(mt, 'lb', pos{j}, -d(i,j)); else, mt = setParam(mt, 'eq', pos{j}, -d(i,j)); end
    end
    mt = setParam(mt, 'obj', 'r_2111', 1);
    s = solveLP(mt, 1);
    other = exch(s.x(exch) > 1e-4 & ~ismember(mt.rxns(exch), products));
    fprintf('point %d (glucose %.2f, NH4 %.3f): growth %.3f (meas %.3f); non-product exports: %s\n', i, d(i,1), d(i,3), ...
        s.x(strcmp(mt.rxns,'r_2111')), d(i,4), strjoin(strcat(mt.rxnNames(other)', '=', arrayfun(@(v) sprintf('%.3f',v), s.x(other)', 'UniformOutput', false)), ', '));
end
end

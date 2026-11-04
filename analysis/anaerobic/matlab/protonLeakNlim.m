function protonLeakNlim()
% protonLeakNlim  Growth benchmark (v9.1.0 + new anaerobic constraints), all
%   32 chemostats: predicted growth with and without the free proton leaks
%   r_1824 (H+[e] -> H+[c]) and r_2129 (H+[c] -> H+[m]), and their pFBA flux.
paths = setupStudy('3.0.0-rc1');
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
base = loadTaggedModel(paths.yeast, '9.1.0');
noLeak = setParam(base, 'ub', {'r_1824','r_2129'}, 0);
figure('Visible','off'); [~, e, g0] = growthNewSim(base); close all
figure('Visible','off'); [~, ~, g1] = growthNewSim(noLeak); close all
o = double(e(:,4));
fprintf('growth R2 with leaks %.4f, without %.4f\n', 1-sum((o-g0(:,4)).^2)/sum((o-mean(o)).^2), 1-sum((o-g1(:,4)).^2)/sum((o-mean(o)).^2));
fprintf('max |growth difference| %.2e\n', max(abs(g0(:,4)-g1(:,4))));
end

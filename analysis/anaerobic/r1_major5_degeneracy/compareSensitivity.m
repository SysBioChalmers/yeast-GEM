function compareSensitivity()
% compareSensitivity  Samples on the optimal face (growth or glucose at
%   0.99999 of the optimum) against samples with 1% slack (0.99), for the 16
%   sensitivity runs (4 variants x Sjoberg, Bakker wild type, Jouhten 2008 and
%   Wasylenko 2014). Writes results/analysis_sensitivity/{exchanges,fluxomics,claims}.tsv.

paths = setupStudy('3.0.0-rc1');
in  = fullfile(paths.results, 'sampling');
out = fullfile(paths.results, 'analysis_sensitivity');
if ~isfolder(out), mkdir(out); end
variants = {'vdev_new','vdev_old','v902_new','v902_old'};
opts = [0.99999 0.99];
load1 = @(v, s, o) load(fullfile(in, sprintf('%s__%s__opt%g.mat', v, s, o)));
q = @(x) prctile(x, [2.5 50 97.5], 2);

ex = {'r_4041','growth'; 'r_1808','glycerol'; 'r_1761','ethanol'; 'r_1672','CO2'; 'r_1714','glucose'};
T = table();
for v = variants
    for s = {'sjoberg','bakker_wt'}
        for o = opts
            R = load1(v{1}, s{1}, o);
            for k = 1:size(ex,1)
                qq = q(abs(R.W(strcmp(R.rxns, ex{k,1}), :)));
                T = [T; table(v, s, o, ex(k,2), qq(1), qq(2), qq(3), 'VariableNames', ...
                    {'variant','scenario','optFrac','quantity','q2_5','median','q97_5'})]; %#ok<AGROW>
            end
        end
    end
end
writetable(T, fullfile(out, 'exchanges.tsv'), 'FileType','text','Delimiter','\t');

F = readtable(fullfile(paths.yeast, 'data', 'physiology', 'flux_data_anaerobic.tsv'), ...
    'FileType','text','ReadVariableNames',false,'Delimiter','\t');
P = table();
for v = variants
    for s = {'Jouhten2008','Wasylenko2014'}
        rows = find(strcmp(F{:,6}, s{1}));
        for o = opts
            R = load1(v{1}, ['flux_' s{1}], o);
            glc = find(strcmp(R.rxns, 'r_1714'));
            c = []; cov = []; w = [];
            for r = rows'
                i = find(strcmp(R.rxns, F{r,7}));
                if isempty(i) || i == glc, continue; end
                x = abs(100 * R.W(i,:) ./ R.W(glc,:)); y = abs(F{r,5}); qq = q(x);
                c(end+1) = crps(x, y); cov(end+1) = y >= qq(1) & y <= qq(3); w(end+1) = qq(3) - qq(1); %#ok<AGROW>
            end
            P = [P; table(v, s, o, numel(c), mean(c), mean(cov), median(w), 'VariableNames', ...
                {'variant','dataset','optFrac','n','crps_mean','coverage95','median_width95'})]; %#ok<AGROW>
        end
    end
end
writetable(P, fullfile(out, 'fluxomics.tsv'), 'FileType','text','Delimiter','\t');

claims = {'r_0713','Mdh1'; 'r_2131','Idp1'; 'r_0658','Idh1/2'; 'r_0471','Gdh1/3'; 'r_0476','Gln1'; ...
          'r_0226','mito ATP synthase'; 'r_0227','PM H+-ATPase'; 'r_0491','Gpd1/2'; 'r_1810','glycine export'};
C = table();
for o = opts
    R = load1('vdev_new', 'sjoberg', o);
    for k = 1:size(claims,1)
        qq = q(R.W(strcmp(R.rxns, claims{k,1}), :));
        C = [C; table(o, claims(k,2), qq(1), qq(2), qq(3), 'VariableNames', {'optFrac','reaction','q2_5','median','q97_5'})]; %#ok<AGROW>
    end
end
writetable(C, fullfile(out, 'claims.tsv'), 'FileType','text','Delimiter','\t');
disp(T); disp(P); disp(C);
end

function c = crps(x, y)
% Empirical continuous ranked probability score, E|X-y| - E|X-X'|/2.
x = sort(x(:)); n = numel(x);
c = mean(abs(x - y)) - sum((2*(1:n)' - n - 1) .* x) / n^2;
end

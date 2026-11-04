function analyzeSampling(optFrac, variants, tag)
% analyzeSampling  Turn the sampling results into the tables for the revision.
%   analyzeSampling(0.99999)
%   analyzeSampling(0.99999, {'v911_new','v911_old','v911c_new','v902_new','v902_old'}, 'v911')
%   variants default to the four v9.1.0/v9.0.2 variants; tag, if given, is
%   appended to the output folder name.
%
%   Reads results/sampling/<variant>__<scenario>__opt<optFrac>.mat (cycle-free
%   samples W) and writes to results/analysis_opt<optFrac>/:
%     fig1_sjoberg.tsv      growth and product exchanges: sampled 2.5/50/97.5%
%                           and pFBA, against the measured value (Figure 1)
%     fluxomics_points.tsv  per measured flux ratio (100*|v/v_glc|): sampled
%                           quantiles, pFBA, CRPS, and whether the measurement
%                           lies in the 95% interval (Figure 3C-D)
%     fluxomics_summary.tsv per variant: mean CRPS, coverage, median interval
%                           width, pFBA MRE; with and without the glucose points
%     knockouts.tsv         Bakker 2000 and Boles 1998 scenarios against data
%     claims.tsv            reactions discussed in the text: sampled interval and
%                           probability of being active in each direction
%     protons.tsv           every reaction moving H+ between compartments

if nargin < 1, optFrac = 0.99999; end
paths = setupStudy('3.0.0-rc1');
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
in  = fullfile(paths.results, 'sampling');
if nargin < 2 || isempty(variants), variants = {'v910_new','v910_old','v902_new','v902_old'}; end
if nargin < 3, tag = ''; end
out = fullfile(paths.results, sprintf('analysis_opt%g', optFrac));
if ~isempty(tag), out = [out '_' tag]; end
if ~isfolder(out), mkdir(out); end
load1 = @(v, s) load(fullfile(in, sprintf('%s__%s__opt%g.mat', v, s, optFrac)));
q = @(x) prctile(x, [2.5 50 97.5], 2);

%% Figure 1: Sjoberg batch
meas = {'r_4041','growth',0.36,0.02; 'r_1808','glycerol',4.5,0.4; 'r_1761','ethanol',31,2; 'r_1672','CO2',38,10};
T = table();
for v = variants
    R = load1(v{1}, 'sjoberg');
    for k = 1:size(meas,1)
        i = strcmp(R.rxns, meas{k,1}); x = abs(R.W(i,:)); qq = q(x);
        T = [T; table(v, meas(k,2), meas{k,3}, meas{k,4}, abs(R.spec.pfba(i)), qq(1), qq(2), qq(3), ...
            meas{k,3} >= qq(1)-meas{k,4} & meas{k,3} <= qq(3)+meas{k,4}, ...
            'VariableNames', {'variant','quantity','measured','sd','pfba','q2_5','median','q97_5','overlaps_within_sd'})]; %#ok<AGROW>
    end
end
writetable(T, fullfile(out, 'fig1_sjoberg.tsv'), 'FileType','text','Delimiter','\t');

%% Figure 3C-D: fluxomics
F = readtable(fullfile(paths.yeast, 'data', 'physiology', 'flux_data_anaerobic.tsv'), ...
    'FileType','text','ReadVariableNames',false,'Delimiter','\t');
sets = unique(F{:,6});
P = table();
for v = variants
    for s = sets'
        R = load1(v{1}, ['flux_' s{1}]);
        rows = find(strcmp(F{:,6}, s{1}));
        glc = find(strcmp(R.rxns, 'r_1714'));
        for r = rows'
            i = find(strcmp(R.rxns, F{r,7}));
            if isempty(i), continue; end            % as in the v9.1.0 benchmark
            x = abs(100 * R.W(i,:) ./ R.W(glc,:));
            y = abs(F{r,5});
            p = abs(100 * R.spec.pfba(i) / R.spec.pfba(glc));
            qq = q(x);
            P = [P; table(v, s, F{r,7}, y, p, qq(1), qq(2), qq(3), crps(x, y), ...
                y >= qq(1) & y <= qq(3), qq(3)-qq(1), strcmp(F{r,7}, 'r_1714'), ...
                'VariableNames', {'variant','dataset','rxn','measured','pfba','q2_5','median','q97_5','crps','in95','width95','is_glucose'})]; %#ok<AGROW>
        end
    end
end
writetable(P, fullfile(out, 'fluxomics_points.tsv'), 'FileType','text','Delimiter','\t');
Sm = table();
for v = variants
    for incl = [true false]
        sel = strcmp(P.variant, v{1}) & (incl | ~P.is_glucose);
        re = abs(P.pfba(sel) - P.measured(sel)) ./ P.measured(sel);
        Sm = [Sm; table(v, incl, nnz(sel), mean(P.crps(sel)), median(P.crps(sel)), mean(P.in95(sel)), ...
            median(P.width95(sel)), mean(re(isfinite(re))), ...
            'VariableNames', {'variant','with_glucose_points','n','crps_mean','crps_median','coverage95','median_width95','pfba_MRE'})]; %#ok<AGROW>
    end
end
writetable(Sm, fullfile(out, 'fluxomics_summary.tsv'), 'FileType','text','Delimiter','\t');

%% Knockouts
[kHead, kRows] = readKnockoutData(fullfile(paths.study, 'data', 'knockout_data.tsv'));
prods = {'r_4041','mu'; 'r_1761','q_ethanol'; 'r_1808','q_glycerol'; 'r_1634','q_acetate'; 'r_2056','q_succinate'; 'r_1672','q_CO2'; 'r_1714','q_glucose'};
scen = {'bakker_wt','wild type'; 'bakker_adh3','adh3'; 'bakker_adh3adh4',''; 'boles_wt',''; 'boles_mae1',''};
Ko = table();
for v = variants
    for s = 1:size(scen,1)
        R = load1(v{1}, scen{s,1});
        for k = 1:size(prods,1)
            m = NaN; sd = NaN;
            c = find(strcmp(kHead, prods{k,2}));
            if ~isempty(scen{s,2}) && ~isempty(c)
                row = find(strcmp(kRows(:,1),'Bakker2000') & strcmp(kRows(:,4),'anaerobic') & strcmp(kRows(:,3), scen{s,2}));
                if ~isempty(row), m = str2double(kRows{row, c}); sd = str2double(kRows{row, c+1}); end
            end
            if ~R.spec.feasible
                Ko = [Ko; table(v, scen(s,1), prods(k,2), m, sd, 0, NaN, NaN, NaN, false, ...
                    'VariableNames', {'variant','scenario','quantity','measured','measured_sd','pfba','q2_5','median','q97_5','feasible'})]; %#ok<AGROW>
                continue
            end
            i = strcmp(R.rxns, prods{k,1}); x = abs(R.W(i,:)); qq = q(x);
            Ko = [Ko; table(v, scen(s,1), prods(k,2), m, sd, abs(R.spec.pfba(i)), qq(1), qq(2), qq(3), true, ...
                'VariableNames', {'variant','scenario','quantity','measured','measured_sd','pfba','q2_5','median','q97_5','feasible'})]; %#ok<AGROW>
        end
    end
end
writetable(Ko, fullfile(out, 'knockouts.tsv'), 'FileType','text','Delimiter','\t');

%% Claims in the text and proton transport (Sjoberg scenario)
claims = {'r_0713','Mdh1 (mitochondrial malate dehydrogenase)'; 'r_0714','Mdh2 (cytosolic, blocked)'; ...
          'r_2131','Idp1 (mitochondrial NADP isocitrate dehydrogenase)'; 'r_0659','Idp2 (cytosolic NADP isocitrate dehydrogenase, blocked)'; ...
          'r_0658','Idh1/2 (mitochondrial NAD isocitrate dehydrogenase)'; ...
          'r_0471','Gdh1/3 (NADP glutamate dehydrogenase)'; 'r_0470','Gdh2 (NAD glutamate dehydrogenase)'; ...
          'r_0476','Gln1 (glutamine synthetase)'; 'r_0472','Glt1 (glutamate synthase)'; ...
          'r_0501','glycine cleavage'; 'r_0487','Gcy1 (glycerol dehydrogenase, NADP)'; ...
          'r_0226','mitochondrial ATP synthase'; 'r_0227','plasma-membrane H+-ATPase'; ...
          'r_1654','ammonium exchange'; 'r_1115','ammonium permease transport'; 'r_0491','Gpd1/2 (glycerol-3-phosphate dehydrogenase)'; ...
          'r_0252','carnitine acetyltransferase'};
C = table(); H = table();
for v = variants
    R = load1(v{1}, 'sjoberg');
    for k = 1:size(claims,1)
        i = find(strcmp(R.rxns, claims{k,1}));
        if isempty(i), continue; end
        x = R.W(i,:); qq = q(x);
        C = [C; table(v, claims(k,1), claims(k,2), R.spec.pfba(i), qq(1), qq(2), qq(3), ...
            mean(x > 1e-6), mean(x < -1e-6), 'VariableNames', ...
            {'variant','rxn','description','pfba','q2_5','median','q97_5','p_forward','p_reverse'})]; %#ok<AGROW>
    end
    % Reactions that move H+ between compartments.
    m = modelFor(paths, v{1});
    hIdx = find(strcmp(m.metNames, 'H+'));
    for i = find(sum(m.S(hIdx,:) ~= 0, 1) >= 2)
        if all(sign(m.S(hIdx, i)) >= 0) || all(sign(m.S(hIdx, i)) <= 0), continue; end
        j = find(strcmp(R.rxns, m.rxns{i}));
        x = R.W(j,:); qq = q(x);
        others = nnz(m.S(:,i)) - nnz(m.S(hIdx,i));
        H = [H; table(v, m.rxns(i), m.rxnNames(i), constructEquations(m, m.rxns(i)), others == 0, ...
            R.spec.pfba(j), qq(1), qq(2), qq(3), mean(abs(x) > 1e-6), 'VariableNames', ...
            {'variant','rxn','name','equation','proton_only','pfba','q2_5','median','q97_5','p_active'})]; %#ok<AGROW>
    end
end
writetable(C, fullfile(out, 'claims.tsv'), 'FileType','text','Delimiter','\t');
writetable(H, fullfile(out, 'protons.tsv'), 'FileType','text','Delimiter','\t');
disp(T); disp(Sm);
end

function [head, rows] = readKnockoutData(file)
% Tab-separated, '#' comment lines; returns the header and the rows as text.
lines = strsplit(fileread(file), {'\n'});
lines = lines(~startsWith(lines, '#') & strlength(strtrim(lines)) > 0);
head = strsplit(strtrim(lines{1}), '\t', 'CollapseDelimiters', false);
rows = repmat({''}, numel(lines)-1, numel(head));
for i = 2:numel(lines)
    f = strsplit(strtrim(lines{i}), '\t', 'CollapseDelimiters', false);
    rows(i-1, 1:numel(f)) = f;
end
end

function c = crps(x, y)
% Empirical continuous ranked probability score, E|X-y| - E|X-X'|/2.
x = sort(x(:)); n = numel(x);
c = mean(abs(x - y)) - sum((2*(1:n)' - n - 1) .* x) / n^2;
end

function m = modelFor(paths, variant)
persistent cache
if isempty(cache), cache = containers.Map(); end
if ~isKey(cache, variant)
    cache(variant) = variantModel(paths, variant);
end
m = cache(variant);
end

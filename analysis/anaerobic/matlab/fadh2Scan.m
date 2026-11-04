function fadh2Scan(outName)
% fadh2Scan  Yeast9.1.0 with the updated anaerobic constraints, with the FADH2
%   coefficient of the anaerobic biomass reaction (0.08 as released) set to a
%   range of values: 13C flux benchmark (MRE without the glucose points),
%   predicted succinate against the measured values in the flux data, Sjoberg
%   glycerol and succinate, and anaerobic frd1 osm1 viability.
%   Writes results/<outName>.tsv (scan) and results/<outName>_succinate.tsv.
if nargin < 1, outName = 'fadh2_scan'; end
paths = setupStudy('3.0.0-rc1');
cd(fullfile(paths.yeast,'code','otherChanges'));
base = variantModel(paths,'v910_new');
bio = strcmp(base.rxns,'r_4041');
[~,mi] = ismember({'s_0689','s_0687','s_0794'}, base.mets);   % FADH2[c], FAD[c], H+[c]
assert(abs(base.S(mi(1),bio) - 0.08) < 1e-9, 'FADH2 coefficient is not 0.08');
vals = [0 0.005 0.01 0.019 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10 0.12 0.16 0.20 0.30];
fid = fopen(fullfile(paths.results,[outName '.tsv']),'w');
fprintf(fid,'fadh2\tflux_MRE_noglc\tflux_MRE_all\tflux_R2\tsuccinate_MRE\tsjoberg_growth\tsjoberg_glycerol\tsjoberg_succinate\tfrd1osm1_rel_growth\n');
fs = fopen(fullfile(paths.results,[outName '_succinate.tsv']),'w');
fprintf(fs,'fadh2\tdataset\tmeasured\tpredicted\n');
for c = vals
    m = base; d = c - 0.08;
    m.S(mi,bio) = m.S(mi,bio) + [d; -d; -2*d];
    figure('Visible','off');
    [R2,mreAll,D,S,N,T] = fluxPredSim(m);
    close all
    keep = ~strcmp(N,'r_1714');
    mre = mean(abs((S(keep)-D(keep))./D(keep)),'omitnan');
    su = strcmp(N,'r_2056');
    sre = mean(abs((S(su)-D(su))./D(su)),'omitnan');
    for k = find(su)
        fprintf(fs,'%g\t%s\t%.4f\t%.4f\n', c, T{k}, D(k), S(k));
    end
    sim = sjobergSim(m);
    s23 = solveLP(setParam(m,'eq','r_1714',-23),1);
    succ = s23.x(strcmp(m.rxns,'r_2056'));
    a20 = setParam(setParam(m,'eq','r_1714',-20),'obj','r_4041',1);
    s0 = solveLP(a20,1); g0 = s0.x(bio);
    k2 = removeGenes(a20,{'YEL047C','YJR051W'},false,false,false); sk = solveLP(k2,1);
    g2 = 0; if ~isempty(sk.x), g2 = sk.x(strcmp(k2.rxns,'r_4041')); end
    fprintf(fid,'%g\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.3f\n', c, mre, mreAll, R2, sre, sim(4), sim(1), succ, g2/g0);
    fprintf('FADH2 %.3f: flux MRE %.4f (all %.4f), succinate MRE %.3f, Sjoberg glycerol %.3f succinate %.4f, frd1 osm1 %.2f\n', ...
        c, mre, mreAll, sre, sim(1), succ, g2/g0);
end
fclose(fid); fclose(fs);
end

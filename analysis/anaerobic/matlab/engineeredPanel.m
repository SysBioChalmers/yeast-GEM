function engineeredPanel(variants)
% engineeredPanel  Anaerobic predictions for engineered strains from five
%   Delft/DTU studies (data/engineered_strain_data.tsv), per model variant:
%   A de Kok 2011: maltose-limited chemostat D 0.05; hydrolysis (maltose-H+
%     symport + maltase) vs phosphorolysis (maltase removed; maltose
%     phosphorylase: maltose + Pi -> glucose + glucose 1-phosphate; the
%     native phosphoglucomutase takes glucose 1-phosphate to glucose 6-P).
%     Biomass yield in g per g glucose equivalent.
%   B Basso 2011: sucrose-limited chemostat D 0.10; extracellular invertase
%     (SUC2) vs intracellular hydrolysis (invertase removed; sucrose-H+
%     symport added). Biomass yield in g per mol sucrose.
%   C Kozak 2014: glucose-limited chemostat D 0.10; wild type vs acs1 acs2 +
%     pyruvate formate lyase (pyruvate + CoA -> acetyl-CoA + formate).
%   D Nissen 2000 (Metab Eng): anaerobic batch; gdh1 (GDH3 also off, repressed
%     on glucose), + GLT1 GLN1 (glutamate synthase opened), + GDH2 (already
%     active). C-mol yields.
%   E Bro 2006: anaerobic batch; + non-phosphorylating NADP-glyceraldehyde-
%     3-phosphate dehydrogenase (GAPN). Yields in g/g glucose.
%   Writes results/engineered_panel.tsv.

if nargin < 1, variants = {'v910_new','v910_old','v902_new','v902_old'}; end
paths = setupStudy('3.0.0-rc1');
cd(fullfile(paths.yeast, 'code', 'otherChanges'));
T = table();
for v = variants
    an0 = variantModel(paths, v{1});
    an0 = setParam(an0, 'obj', 'r_4041', 1);
    g1p = an0.metNames{find(contains(an0.metNames, 'glucose 1-phosphate') & strcmp(an0.comps(an0.metComps), 'c'), 1)};

    % A: maltose chemostat
    base = setParam(setParam(an0, 'lb', 'r_1714', 0), 'lb', 'r_1931', -1000);
    hyd = base;
    pho = addEq(setParam(base, 'eq', 'r_0198', 0), 'MP', 'maltose phosphorylase', ...
        ['maltose[c] + phosphate[c] => D-glucose[c] + ' g1p '[c]']);
    for c = {'hydrolytic (IMZ199)', hyd; 'phosphorolytic (IMZ226)', pho}'
        q = chemostat(c{2}, 0.05, 'r_1931');
        T = addRow(T, v{1}, 'de Kok 2011', c{1}, 'Y_biomass (g/g glucose eq)', 0.05 / (q * 2 * 0.18016));
    end
    q = chemostat(an0, 0.05, 'r_1714');
    T = addRow(T, v{1}, 'de Kok 2011', 'glucose reference', 'Y_biomass (g/g glucose eq)', 0.05 / (q * 0.18016));

    % B: sucrose chemostat
    base = setParam(setParam(an0, 'lb', 'r_1714', 0), 'lb', 'r_2058', -1000);
    ext = base;
    int = addEq(setParam(base, 'eq', 'r_1024', 0), 'SUCt', 'sucrose-H+ symport', ...
        'sucrose[e] + H+[e] => sucrose[c] + H+[c]');
    for c = {'SUC2 (CEN.PK113-7D)', ext; 'iSUC2 (IMI056)', int}'
        q = chemostat(c{2}, 0.10, 'r_2058');
        T = addRow(T, v{1}, 'Basso 2011', c{1}, 'Y_biomass (g/mol sucrose)', 0.10 / q * 1000);
    end

    % C: acs1 acs2 + PFL, glucose chemostat
    pfl = addEq(removeGenes(an0, {'YAL054C','YLR153C'}), 'PFL', 'pyruvate formate lyase', ...
        'pyruvate[c] + coenzyme A[c] => acetyl-CoA[c] + formate[c]');
    for c = {'CEN.PK113-7D', an0; 'acs1 acs2 + PFL (IMZ383)', pfl}'
        [q, s, m] = chemostat(c{2}, 0.10, 'r_1714');
        T = addRow(T, v{1}, 'Kozak 2014', c{1}, 'Y_biomass (g/g glucose)', 0.10 / (q * 0.18016));
        if ~isempty(s), T = addRow(T, v{1}, 'Kozak 2014', c{1}, 'q_formate', s.x(strcmp(m.rxns, 'r_1793'))); end
    end

    % D: ammonium assimilation, batch, C-mol yields
    gdh = removeGenes(an0, {'YOR375C','YAL062W'});
    for c = {'TN1 wild type', an0; 'TN9 gdh1', gdh; 'TN19 gdh1 +GLT1 +GLN1', setParam(gdh, 'ub', 'r_0472', 1000)}'
        s = batch(c{2});
        cm = @(r, nC) nC * s.x(strcmp(c{2}.rxns, r)) / (6 * 20);
        if isempty(s.x), vals = [0 NaN NaN NaN]; else
            mu = s.x(strcmp(c{2}.rxns, 'r_4041'));   % biomass taken as 24.6 g per C-mol
            vals = [mu, mu / 0.0246 / (6 * 20), cm('r_1761', 2), cm('r_1808', 3)]; end
        T = addRow(T, v{1}, 'Nissen 2000 ME', c{1}, 'mu (glucose 20)', vals(1));
        T = addRow(T, v{1}, 'Nissen 2000 ME', c{1}, 'Y_biomass (Cmol/Cmol)', vals(2));
        T = addRow(T, v{1}, 'Nissen 2000 ME', c{1}, 'Y_ethanol (Cmol/Cmol)', vals(3));
        T = addRow(T, v{1}, 'Nissen 2000 ME', c{1}, 'Y_glycerol (Cmol/Cmol)', vals(4));
    end

    % E: GAPN, batch
    gapn = addEq(an0, 'GAPN', 'non-phosphorylating glyceraldehyde-3-phosphate dehydrogenase (NADP)', ...
        'glyceraldehyde 3-phosphate[c] + NADP(+)[c] + H2O[c] => 3-phosphonato-D-glycerate(3-)[c] + NADPH[c] + 2 H+[c]');
    for c = {'reference', an0; 'gapN (CPB1001)', gapn}'
        s = batch(c{2});
        y = @(r, mw) s.x(strcmp(c{2}.rxns, r)) * mw / (20 * 180.156);
        T = addRow(T, v{1}, 'Bro 2006', c{1}, 'Y_biomass (g/g)', s.x(strcmp(c{2}.rxns, 'r_4041')) / (20 * 0.180156));
        T = addRow(T, v{1}, 'Bro 2006', c{1}, 'Y_glycerol (g/g)', y('r_1808', 92.09));
        T = addRow(T, v{1}, 'Bro 2006', c{1}, 'Y_ethanol (g/g)', y('r_1761', 46.07));
    end
    fprintf('%s done\n', v{1});
end
writetable(T, fullfile(paths.results, 'engineered_panel.tsv'), 'FileType','text','Delimiter','\t');
disp(T);
end

function m = addEq(m, id, name, eq)
r.rxns = {id}; r.rxnNames = {name}; r.equations = {eq};
r.lb = 0; r.ub = 1000;
m = addRxns(m, r, 3, [], false);
end

function [q, s, m] = chemostat(m, D, exch)
% Growth fixed at D; the least substrate uptake that supports it (pFBA).
m = setParam(m, 'eq', 'r_4041', D);
m = setParam(m, 'lb', exch, -1000);       % substrate uptake limited only by the objective
m = setParam(m, 'obj', exch, 1);          % maximise the (negative) exchange
s = solveLP(m, 1);
if isempty(s.x), q = NaN; s = []; return; end
q = -s.x(strcmp(m.rxns, exch));
end

function s = batch(m)
m = setParam(setParam(m, 'eq', 'r_1714', -20), 'obj', 'r_4041', 1);
s = solveLP(m, 1);
end

function T = addRow(T, v, src, strain, q, val)
T = [T; table({v}, {src}, {strain}, {q}, val, 'VariableNames', {'variant','source','strain','quantity','predicted'})];
end

function [modelAn, base, anaerobicFn, post] = variantModel(paths, variant)
% variantModel  Model for a variant name '<release>_<constraints>'.
%   release      v902, vdev, ... (yeast-GEM tag), or vdev (the checkout: yeast-GEM 9.1.2, develop with #423 and #426 merged)
%   constraints  'new' (applyEnvironment(model,'anaerobic')) or 'old' (the
%                constraints before 9.1.0, anaerobicModel of this checkout)
%   Returns the anaerobic model, the aerobic base model, the anaerobic function and the step applied after it; the
%   growth benchmark takes base, anaerobicFn is implied by the benchmark copy
%   and post is passed to it.
%   Needs the current folder to be code/otherChanges (anaerobicModel).
%   Optional '+<modifier>' suffixes add a step after the anaerobic function
%   (see modifierStep), e.g. 'vdev_new+mdh1old'.

parts = strsplit(variant, '+');
tok = regexp(parts{1}, '^v(dev|\d\d\d)_(new|old)$', 'tokens', 'once');
assert(~isempty(tok), 'Unknown variant %s', variant);
base = loadTaggedModel(paths.yeast, regexprep(tok{1}, '^(\d)(\d)(\d)$', '$1.$2.$3'));
post = @(m) m;
if strcmp(tok{1}, 'dev') && strcmp(tok{2}, 'old')
    post = @openThiamine;   % the curated model needs thiamine uptake under anaerobic conditions
end
for k = 2:numel(parts)
    post = @(m) modifierStep(post(m), parts{k});
end
if strcmp(tok{2}, 'new'), anaerobicFn = @(m) applyEnvironment(m, 'anaerobic'); else, anaerobicFn = @anaerobicModel; end
modelAn = post(anaerobicFn(base));
end

function model = modifierStep(model, name)
switch name
    case 'mdh1old'   % restriction of the original constraints: mitochondrial Mdh1 malate -> OAA only
        model.lb(strcmp(model.rxns, 'r_0713')) = 0;
    case {'fadh2x019','fadh2x0'}   % FADH2 production in biomass: 0.08 -> 0.019 (Cys ceiling) or 0
        target = 0.019*strcmp(name,'fadh2x019');
        bio = strcmp(model.rxns,'r_4041');
        [~,mi] = ismember({'s_0689','s_0687','s_0794'}, model.mets);   % FADH2[c], FAD[c], H+[c]
        assert(abs(model.S(mi(1),bio) - 0.08) < 1e-9, 'FADH2 coefficient is not 0.08');
        d = target - 0.08;
        model.S(mi,bio) = model.S(mi,bio) + [d; -d; -2*d];
    case 'glt1open'   % glutamate synthase unblocked
        model.ub(strcmp(model.rxns,'r_0472')) = 1000;
    case 'glt1cap'    % Glt1 capacity from Table 1: 0.75 mg/gprot x 42 U/mg x 0.46 gprot/gDW = 0.87 mmol/gDW/h
        model.ub(strcmp(model.rxns,'r_0472')) = 0.87;
    case 'gcvopen'
        model.ub(ismember(model.rxns,{'r_0501','r_0507','r_0509'})) = 1000;
    case 'gcvcap'     % GCV capacity from Table 1: 0.14 mg/gprot x 1.25 U/mg x 0.46 gprot/gDW = 0.005 mmol/gDW/h
        model.ub(ismember(model.rxns,{'r_0501','r_0507','r_0509'})) = 0.005;
    case 'gcy1rev'    % Gcy1 reversible (9.1.0: DHA + NADPH -> glycerol only)
        model.lb(strcmp(model.rxns,'r_0487')) = -1000;
    case 'gcy1revcap' % Gcy1 reversible, capacity from Table 1: 0.32 mg/gprot x 12 U/mg x 0.46 gprot/gDW = 0.106 mmol/gDW/h
        model.lb(strcmp(model.rxns,'r_0487')) = -0.106;
        model.ub(strcmp(model.rxns,'r_0487')) = min(model.ub(strcmp(model.rxns,'r_0487')), 0.106);
    case 'undoTetraenol'   % updated function blocks ergosta-5,7,22,24(28)-tetraen-3beta-ol uptake
        model.lb(strcmp(model.rxns,'r_2137')) = -1000;
    case 'undoMdh2'        % updated function blocks Mdh2
        model.lb(strcmp(model.rxns,'r_0714')) = -1000; model.ub(strcmp(model.rxns,'r_0714')) = 1000;
    case 'undoIdp2'        % updated function blocks Idp2
        model.lb(strcmp(model.rxns,'r_0659')) = -1000; model.ub(strcmp(model.rxns,'r_0659')) = 1000;
    case 'undoAA'          % aerobic instead of anaerobic amino-acid ratio
        model = changeAminoAcidRatio(model,true);
    case 'oldGAM'          % original function: GAM 30.49, NGAM 0, protein 0.461 (Nissen 1997)
        model = changeGAM(model,30.49,0);
        model = scaleBioMass(model,'protein',0.461,'carbohydrate',false);
    case 'oldCofactor'     % original function: NAD, NADH, NADP, NADPH, CoA removed from the cofactor pseudoreaction
        mi = ismember(model.mets,{'s_1198','s_1203','s_1207','s_1212','s_0529'});
        model.S(mi,strcmp(model.rxns,'r_4598')) = 0;
    case 'gcy1ox'     % Gcy1 glycerol -> DHA only (direction documented in vivo)
        model.lb(strcmp(model.rxns,'r_0487')) = -1000;
        model.ub(strcmp(model.rxns,'r_0487')) = 0;
    otherwise
        error('Unknown modifier %s', name);
end
end

function model = openThiamine(model)
% Thiamine uptake, as in data/conditions/anaerobic.yml: the original constraints
% do not provide it, and the oxygen-dependent HMP-P synthase is blocked anaerobically.
model.lb(strcmp(model.rxns, 'r_2067')) = -1000;
end

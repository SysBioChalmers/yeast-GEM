function biomassDR()
% biomassDR  Degree of reduction of the biomass made by the model, per C-mol,
%   from the exchange fluxes (all mass enters and leaves through exchanges):
%   DR_biomass = -sum(v_ex .* DR_met) / -sum(v_ex .* C_met), with DR of C 4, H 1,
%   O -2, N -3 (NH3 reference), S 6, P 5 and -1 per positive charge (Roels 1983).
%   Aerobic minimal medium at growth 0.1 (glucose minimised), with and without
%   the 3 NADPH/gDW that v9_1_0.m added to the biomass reaction to reach the
%   4.2 per C-mol of Lange and Heijnen (2001).
paths = setupStudy('3.0.0-rc1');
cd(fullfile(paths.yeast,'code','otherChanges'));
base = loadTaggedModel(paths.yeast,'9.1.0');
noDR = base;
bio = strcmp(noDR.rxns,'r_4041');
[~,o] = ismember({'s_1212','s_1207','s_0794'}, noDR.mets);
noDR.S(o,bio) = noDR.S(o,bio) - [-3; 3; -3];
b902 = loadTaggedModel(paths.yeast,'9.0.2');
cases = {{'9.0.2', b902, false}, {'9.1.0', base, false}, {'9.1.0 without biomass NADPH', noDR, false}, ...
         {'9.1.0 anaerobic (updated constraints)', base, true}, {'9.1.0 anaerobic without biomass NADPH', noDR, true}};
for c = cases
    m = c{1}{2};
    m = minimal_Y6_local(m);
    if c{1}{3}, m = anaerobicModel(m); end
    m = setParam(m,'lb','r_1714',-10);
    m = setParam(m,'lb','r_4041',0.1);
    m = setParam(m,'obj','r_1714',1);       % minimise glucose uptake (maximise the negative flux)
    s = solveLP(m,1);
    ex = find(sum(m.S~=0,1)==1)';            % exchange reactions
    ex = ex(~strcmp(m.rxns(ex),'r_2111'));   % not the biomass drain itself
    [C,DR] = deal(zeros(numel(ex),1));
    for k = 1:numel(ex)
        met = find(m.S(:,ex(k)));
        [C(k),DR(k)] = elementDR(m.metFormulas{met}, m.metCharges(met));
        C(k) = C(k)*m.S(met,ex(k)); DR(k) = DR(k)*m.S(met,ex(k));
    end
    v = s.x(ex);
    cmol = -sum(v.*C); dr = -sum(v.*DR);
    fprintf('%-32s biomass: %.3f C-mol/gDW, degree of reduction %.3f per C-mol\n', c{1}{1}, cmol/0.1/1000, dr/cmol);
end
end

function m = minimal_Y6_local(m)
here = cd(fullfile(fileparts(which('changeGAM')),'..','modelCuration'));
back = onCleanup(@() cd(here));
m = minimal_Y6(m);
end

function [C,DR] = elementDR(f, charge)
w = struct('C',4,'H',1,'O',-2,'N',-3,'S',6,'P',5);
tok = regexp(f,'([A-Z][a-z]?)(\d*)','tokens');
C = 0; DR = 0;
for t = tok
    n = 1; if ~isempty(t{1}{2}), n = str2double(t{1}{2}); end
    if strcmp(t{1}{1},'C'), C = C + n; end
    if isfield(w,t{1}{1}), DR = DR + n*w.(t{1}{1}); end
end
if ~isnan(charge), DR = DR - charge; end
end

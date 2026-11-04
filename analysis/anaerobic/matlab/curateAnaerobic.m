function model = curateAnaerobic(model, stage)
% curateAnaerobic  Proposed curation package (not in any yeast-GEM release).
%   Only reaction directions and presence change; no flux values are set.
%   stage 'base'       applied to the model before the anaerobic function:
%     - Lot6 (YLR011W) NAD(P)H:FMN reductase on free FMN (r_0441, r_0442)
%       removed: with the FMNH2-linked fumarate reductase it forms a net
%       NADH:fumarate oxidoreductase not seen in vivo; Frd1/Osm1 re-oxidise
%       FADH2 (Camarasa 2007, Kim 2018).
%     - Iah1 (YOR126C) ethyl-ester esterases hydrolysis only (Fukuda 1998,
%       2000; Lilly 2006); the reactions are written as ester synthesis.
%   stage 'anaerobic'  applied after the anaerobic function:
%     - Thi5-family HMP-P synthase (r_1603) needs O2 (Coquille 2012, Mondal
%       2023): blocked; thiamine uptake opened (Verduyn medium contains it).

switch stage
    case 'base'
        model = setParam(model, 'ub', {'r_0441','r_0442'}, 0);
        model = setParam(model, 'ub', {'r_4711','r_4712','r_4715','r_4716','r_4717','r_4720','r_4721'}, 0);
    case 'anaerobic'
        model = setParam(model, 'ub', 'r_1603', 0);
        model = setParam(model, 'lb', 'r_2067', -1000);
    otherwise
        error('Unknown stage %s', stage);
end
end

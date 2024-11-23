function model = anaerobicModel(model,recalBiomassMW)
% anaerobicModel
%   Constrains yeast-GEM to anaerobic conditions. By default yeast-GEM aims
%   to represent aerobic metabolism (particulary with glucose as carbon
%   source). Here, various exchange reactions and a few selected
%   intracellular reactions are enabled/disabled to yield a model that is
%   able to reach similar exchange rates as measured.
%
%   This function was updated as part of release v9.1.0.
%
% Input:
%   model           model structure, which is aerobic by default
%   recalBiomassMW  logical, should the biomass molecular weight be
%                   recalculated before fitting it to 1 g/gDCW. This option
%                   requires COBRA Toolbox to be installed. (optional,
%                   default false, a precalculated value will instead be
%                   used)
%
% Output:
%   model           model structure, modified to match anaerobic conditions
%
% Usage: model = anaerobicModel(model)

if nargin<2
    recalBiomassMW = false;
end

GAM   = 55.2; % Remains unchanged, line can be removed
P     = 0.46;  %Data from Nissen et al. 1997
NGAM  = 1; % 0.7 in aerobic, should it not remain unchanged?

model = changeGAM(model,GAM,NGAM);

%% Set environmental conditions
% Remove heme a from the cofactor pseudoreaction (part of biomass)
hemeIdx = getIndexes(model,'s_3714','mets');
cofacIdx = getIndexes(model,'r_4598','rxns');
model.S(hemeIdx,cofacIdx) = 0;

% Change exchange reactions (block O2 uptake and allow sterol and fatty
% acid exchanges, as these are essential supplements for anaerobic growth).
model.lb(strcmp(model.rxns,'r_1992')) = 0;        %O2
model.lb(strcmp(model.rxns,'r_1757')) = -1000;    %ergosterol
model.lb(strcmp(model.rxns,'r_1915')) = -1000;    %lanosterol
model.lb(strcmp(model.rxns,'r_2106')) = -1000;    %zymosterol
model.lb(strcmp(model.rxns,'r_2134')) = -1000;    %14-demethyllanosterol
model.lb(strcmp(model.rxns,'r_1994')) = -1000;    %palmitoleate
model.lb(strcmp(model.rxns,'r_2189')) = -1000;    %oleate
% Enable uptake of vitamins for NAD(P)H and CoA synthesis
model.lb(strcmp(model.rxns,'r_1967')) = -1000;    %nicotinate
model.lb(strcmp(model.rxns,'r_1548')) = -1000;    %(R)-pantothenate

%% Degree of reduction of biomass
% To align the degree of reduction of S. cerevisiae biomass to the
% published value of 4.2 /Cmol (Lange and Heijnen, 2001, 10.1002/bit.10054)

DR = 3; % 3mmol (g CDW)−1s
metIdx = getIndexes(model,{'s_1212','s_1207','s_0794'},'mets'); % NADPH[c], NADP[c], H+[c]
bioIdx = getIndexes(model,'r_4041','rxns');

currCoeff = full(model.S(metIdx,bioIdx)); % Gather the current coefficients
model.S(metIdx,bioIdx) = currCoeff + [-DR; +DR; +DR];

%% Curations that are required to reach correct metabolic phenotypes during
% anaerobic batch growth on minimal glucose media

% Block MDH2. Involved in growth on two-carbon substrates. Down regulated
% and proteolytically degraded during growth on glucose (Hung et al (2004)
% 10.1074/jbc.M404544200). It is strongly repressed in transcriptome (Tai
% et al (2005) 10.1074 /jbc.M410573200) and not detected in proteome
% (Sjöberg et al (2023) 10.1016/j.ymben.2024.01.007).
model = setParam(model,'eq','r_0714',0);

% Make GCY1 irreversible. Has a positive DeltaGo' (+20.9) and is part of a
% transhydrogenase cycle (NADH -> NADPH) at the cost of one ATP. High
% cytosolic NADPH/NADP ratio makes it thermodynamically infeasible that it
% runs in reverse direction. 
model = setParam(model,'ub','r_0487',0);

% Block IDP2.  It is strongly repressed in transcriptome (Tai et al (2005)
% 10.1074 /jbc.M410573200) and not detected in proteome (Sjöberg et al
% (2023) 10.1016/j.ymben.2024.01.007).
model = setParam(model,'eq',{'r_0659'},0);

%% Fumarate reductase is required for to recycle FADH2 dereived from
% disulphide bound formation in anaerobic conditions through Ero1 (Camarasa
% et al (2007) 10.1002/yea.1467; Kim et al (2018) 10.1038/s41467-018-07285-9). 

FADH2_prod=0.08;
metIdx = getIndexes(model,{'s_0689','s_0687','s_0794'},'mets'); % FADH2[c], FAD[c], H+[c]

currCoeff = full(model.S(metIdx,bioIdx)); % Gather the current coefficients
model.S(metIdx,bioIdx) = currCoeff + [FADH2_prod; -FADH2_prod; -2*FADH2_prod];

%% Adjust Mw of biomass to 1000
% This function requires COBRA Toolbox to be installed.
if recalBiomassMW
    Biomass_MW = computeMetFormulae(model,'metMwRange','s_0450','fillMets','none','printLevel',0);
else
    % Default value, calculate for yeast-GEM 9.1.0
    Biomass_MW = [954.343535827233, 954.343535827305];
end
model.S(:,bioIdx) = model.S(:,bioIdx)*1000/mean(Biomass_MW);
model.S(strcmp(model.mets,'s_0450'),bioIdx) = 1;
end

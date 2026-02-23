function model = anaerobicModel(model)
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
%   model           yeast-GEM model structure, which is aerobic by default
%
% Output:
%   model           model structure, modified to match anaerobic conditions
%
% Usage: model = anaerobicModel(model)

%% Set environmental conditions
% Remove heme a from the cofactor pseudoreaction (part of biomass)
hemeIdx = getIndexes(model,'s_3714','mets');
cofacIdx = getIndexes(model,'r_4598','rxns');
model.S(hemeIdx,cofacIdx) = 0;

model = changeAminoAcidRatio(model,2);

% Change exchange reactions (block O2 uptake and allow sterol and fatty
% acid exchanges, as these are essential supplements for anaerobic growth).
model.lb(strcmp(model.rxns,'r_1992')) = 0;        %O2
model.lb(strcmp(model.rxns,'r_1757')) = -1000;    %ergosterol
model.lb(strcmp(model.rxns,'r_1915')) = -1000;    %lanosterol
model.lb(strcmp(model.rxns,'r_2106')) = -1000;    %zymosterol
model.lb(strcmp(model.rxns,'r_2134')) = -1000;    %14-demethyllanosterol
model.lb(strcmp(model.rxns,'r_1994')) = -1000;    %palmitoleate
model.lb(strcmp(model.rxns,'r_2189')) = -1000;    %oleate
% NEW: remove this due to NADH recycling to ergosterol
model.lb(strcmp(model.rxns,'r_2137')) = 0;    %ergosta-5,7,22,24(28)-tetraen-3beta-ol
% Enable uptake of vitamins for NAD(P)H and CoA synthesis
model.lb(strcmp(model.rxns,'r_1967')) = -1000;    %nicotinate
model.lb(strcmp(model.rxns,'r_1548')) = -1000;    %(R)-pantothenate

%% Curations that are required to reach correct metabolic phenotypes during
% anaerobic batch growth on minimal glucose media

% Block MDH2. Involved in growth on two-carbon substrates. Down regulated
% and proteolytically degraded during growth on glucose (Hung et al (2004)
% 10.1074/jbc.M404544200). It is strongly repressed in transcriptome (Tai
% et al (2005) 10.1074/jbc.M410573200) and not detected in proteome
% (Sjöberg et al (2023) 10.1016/j.ymben.2024.01.007).
model = setParam(model,'eq','r_0714',0);
 
% Block IDP2.  It is strongly repressed in transcriptome (Tai et al (2005)
% 10.1074 /jbc.M410573200) and not detected in proteome (Sjöberg et al
% (2023) 10.1016/j.ymben.2024.01.007).
model = setParam(model,'eq',{'r_0659'},0);

%% Fumarate reductase is required to recycle FADH2 derived from disulphide
% bound formation by growth in anaerobic conditions through Ero1 (Camarasa
% et al (2007) 10.1002/yea.1467; Kim et al (2018) 10.1038/s41467-018-07285-9). 

FADH2_prod=0.08;
metIdx = getIndexes(model,{'s_0689','s_0687','s_0794'},'mets'); % FADH2[c], FAD[c], H+[c]
bioIdx = getIndexes(model,'r_4041','rxns');

currCoeff = full(model.S(metIdx,bioIdx)); % Gather the current coefficients
model.S(metIdx,bioIdx) = currCoeff + [FADH2_prod; -FADH2_prod; -2*FADH2_prod];
end

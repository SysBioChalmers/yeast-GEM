%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% anaerobicModel.m
% Converts model to anaerobic
%
% Benjamín J. Sánchez. Last edited: 2017-10-31
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1st change: Changes media to anaerobic (no O2 uptake and allows sterol
%            and fatty acid exchanges)
model.lb(strcmp(model.rxns,'r_1992')) = 0;        %O2
model.lb(strcmp(model.rxns,'r_1757')) = -1000;    %ergosterol
model.lb(strcmp(model.rxns,'r_1915')) = -1000;    %lanosterol
model.lb(strcmp(model.rxns,'r_1994')) = -1000;    %palmitoleate
model.lb(strcmp(model.rxns,'r_2106')) = -1000;    %zymosterol
model.lb(strcmp(model.rxns,'r_2134')) = -1000;    %14-demethyllanosterol
model.lb(strcmp(model.rxns,'r_2137')) = -1000;    %ergosta-5,7,22,24(28)-tetraen-3beta-ol
model.lb(strcmp(model.rxns,'r_2189')) = -1000;    %oleate

%2nd change: Removes the requirement of heme a in the biomass equation
%            (not used under aerobic conditions)
model.S(strcmp(model.mets,'s_3714[c]'),strcmp(model.rxns,'r_4041')) = 0;

%3rd change: Blocked pathways for proper glycerol production
%Block oxaloacetate-malate shuttle (not present in anaerobic conditions)
model.lb(strcmp(model.rxns,'r_0713')) = 0; %Mithocondria
model.lb(strcmp(model.rxns,'r_0714')) = 0; %Cytoplasm
%Block glycerol dehydroginase (only acts in microaerobic conditions)
model.ub(strcmp(model.rxns,'r_0487')) = 0;
%Block 2-oxoglutarate + L-glutamine -> 2 L-glutamate (alternative pathway)
model.ub(strcmp(model.rxns,'r_0472')) = 0;

%4th change: Refit GAM and NGAM to exp. data, change biomass composition
GAM   = 30.49;  %Data from Nissen et al. 1997
P     = 0.461;  %Data from Nissen et al. 1997
NGAM  = 0;      %Refit done in Jouthen et al. 2012
model = changeBiomass(model,P,GAM,NGAM);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modelCorrections.m
% Corrects various issues in yeast7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Correct glucan coefficients in biomass reaction
model.S(strcmp(model.mets,'s_0002[c]'), strcmp(model.rxns,'r_4041')) = 0;
model.S(strcmp(model.mets,'s_0001[ce]'),strcmp(model.rxns,'r_4041')) = -0.8506;
model.S(strcmp(model.mets,'s_0004[ce]'),strcmp(model.rxns,'r_4041')) = -0.2842;

%Correctly represent proton balance inside cell
model.lb(strcmp(model.rxns,'r_1824')) = 0;  %Block free H+ export
model.ub(strcmp(model.rxns,'r_1250')) = 0;  %Block free putrescine export
model.ub(strcmp(model.rxns,'r_1259')) = 0;  %Block free spermidine export

%Changes in Ox.Pho.
%COMPLEX III: H+ pumping corrected for eff P/O ratio:
eff = 0.633;     %Growth in glucose S.cerevisiae Verduyn et al. 1991
model.S(strcmp(model.mets,'s_0799[m]'),strcmp(model.rxns,'r_0439')) = -2*eff;
model.S(strcmp(model.mets,'s_0794[c]'),strcmp(model.rxns,'r_0439')) = +4*eff;
%COMPLEX IV: H+ pumping corrected for eff P/O ratio:
model.S(strcmp(model.mets,'s_0799[m]'),strcmp(model.rxns,'r_0438')) = -8*eff;
model.S(strcmp(model.mets,'s_0794[c]'),strcmp(model.rxns,'r_0438')) = +4*eff;
%COMPLEX IV: Normalize rxn by the number of ferrocytochromes c:
rxn_pos            = strcmp(model.rxns,'r_0438');
ferro_S            = abs(model.S(strcmp(model.mets,'s_0710[m]'),rxn_pos));
model.S(:,rxn_pos) = model.S(:,rxn_pos)./ferro_S;
%COMPLEX V: For 1 ATP 3 H+ are needed, not 4:
model.S(strcmp(model.mets,'s_0799[m]'),strcmp(model.rxns,'r_0226')) = +2;
model.S(strcmp(model.mets,'s_0794[c]'),strcmp(model.rxns,'r_0226')) = -3;

%Refit GAM and NGAM to exp. data, change biomass composition
GAM   = 40.8;   %Data from Van Hoek at al. 1998
P     = 0.4266; %Data from Van Hoek at al. 1998
NGAM  = 0.7;    %Refit done in Nilsson & Nielsen 2016
model = changeBiomass(model,P,GAM,NGAM);

%Remove unused biomass pseudoreactions:
model = removeRxns(model,{'r_2110','r_2133'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
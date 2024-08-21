function model = anaerobicModel(model)
% anaerobicModel
%   Converts model to anaerobic.
%
%   Inputs: model           (struct) aerobic model
%   Output: model           (struct) anaerobic model
%
%   Usage: model = anaerobicModel(model)
%


GAM   = 55.2;


P     = 0.461;  %Data from Nissen et al. 1997
NGAM  = 1;    

model = changeGAM(model,GAM,NGAM);
model = scaleBioMass(model,'protein',P,'carbohydrate',false);

% %2nd change: Removes the requirement of heme a, NAD(PH), coenzyme A in the biomass equation
% %            (not used under anaerobic conditions)

% I disagree. Please do not remove NAD(P)H and CoA, include the exhange
% reactions for the relevant vitamins instead. See below /Gustav
mets = {'s_3714'}; %,'s_1198','s_1203','s_1207','s_1212','s_0529'};

[~,met_index] = ismember(mets,model.mets);
model.S(met_index,strcmp(model.rxns,'r_4598')) = 0;

%3rd change: Changes media to anaerobic (no O2 uptake and allows sterol
%            and fatty acid exchanges)
model.lb(strcmp(model.rxns,'r_1992')) = 0;        %O2
model.lb(strcmp(model.rxns,'r_1757')) = -1000;    %ergosterol
model.lb(strcmp(model.rxns,'r_1915')) = -1000;    %lanosterol
model.lb(strcmp(model.rxns,'r_1994')) = -1000;    %palmitoleate
model.lb(strcmp(model.rxns,'r_2106')) = -1000;    %zymosterol
model.lb(strcmp(model.rxns,'r_2134')) = -1000;    %14-demethyllanosterol
%remove this due to NADH recycling to ergosterol
model.lb(strcmp(model.rxns,'r_2137')) = 0;    %ergosta-5,7,22,24(28)-tetraen-3beta-ol 
model.lb(strcmp(model.rxns,'r_2189')) = -1000;    %oleate



% Added exchange of vitamins enabling NAD(P)H and CoA syntheis in anaerobic
% conditions /Gustav
model.lb(strcmp(model.rxns,'r_1967')) = -1000;    %nicotinate
model.lb(strcmp(model.rxns,'r_1548')) = -1000;    %(R)-pantothenate


%% Changes to the model that give correct phenotype for anaerobic batch growht on minimal glucose media

% Inhibit MDH2 during excess glucose = anaerobic conditions.
model = setParam(model,'lb','r_0714',0);
model = setParam(model,'ub','r_0714',0);

%Glutamate synthae repressed in excess nitrogen
model.ub(strcmp(model.rxns,'r_0472'))=0;

% GCY1 has a positive DeltaG and is part of a transhydrogenase cycle NADH -> NADPH
model.ub(strcmp(model.rxns,'r_0487')) = 0; 

%glycine cleavage respressed by presence of glucose
model.ub(strcmp(model.rxns,'r_0501'))=0; %glycine cleavage, mitochondrion
model.lb(strcmp(model.rxns,'r_0501'))=0;
model.ub(strcmp(model.rxns,'r_0507'))=0; %glycine cleavage complex (lipoylprotein), mitochondrion
model.lb(strcmp(model.rxns,'r_0507'))=0;
model.ub(strcmp(model.rxns,'r_0509'))=0; %glycine cleavage complex (lipoamide), mitochondrion
model.lb(strcmp(model.rxns,'r_0509'))=0;

% MAE1 and IDP are likely not major mitochondrial NADPH sources. 
% model.ub(strcmp(model.rxns,'r_0719')) = 0; % malic enzyme (MAE1), mitochondrion
% model.ub(strcmp(model.rxns,'r_2131')) = 0; % isocitrate dehydrogenase (IDP1), mitochondrion
% %model.ub(strcmp(model.rxns,'r_0659')) = 0; % isocitrate dehydrogenase (IDP2), mitochondrion

%=========================================================================
%Speculative, this cleans up alternate sources of mitochondrial pyruvate
%(not from the transporter). No effect on glycerol production
%model.ub(strcmp(model.rxns,'r_0718')) = 0; % malic enzyme (MAE1), mitochondrion, NADH reaction, acts as major mitochondrial pyruvate source
%model.ub(strcmp(model.rxns,'r_4701')) = 0; % IRC7, Cysteine desulphydrase, enables growth on cysteine as nitrogen source




% heme a[c]	heme a		bigg.metabolite/hemeA;chebi/CHEBI:24479;kegg.compound/C15670;metanetx.chemical/MNXM53309;sbo/SBO:0000247	C49H55FeN4O6		c	s_3714	-3% cofactor[c]	cofactor		sbo/SBO:0000649			c	s_4205	#NUM!
% r_4598	cofactor pseudoreaction	0.00019 coenzyme A[c] + 1e-05 FAD[c] + 0.00265 NAD[c] + 0.00015 NADH[c] + 0.00057 NADP(+)[c] + 0.0027 NADPH[c] + 0.00099 riboflavin[c] + 1.2e-06 TDP[c] + 6.34e-05 THF[c] + 1e-06 heme a[c] => cofactor[c]
% model.S(find(strcmp(model.mets,'s_3714')),strcmp(model.rxns,'r_4598'))=0;

% A reaction converting NADH to NAD + at 3 mmol (g CDW)−1 was coupled to the growth reaction to give the correct ratio of glycerol production to 
% glucose consumption, and to align the degree of reduction (as defined by (Heijnen, 1994)) of the modeled biomass to a published value for 
% degree of reduction of S. cerevisiae biomass at 4.2 C-mol−1 (Lange and Heijnen, 2001).
% NAD[c]	NAD		bigg.metabolite/nad;chebi/CHEBI:57540;kegg.compound/C00003;metanetx.chemical/MNXM8;sbo/SBO:0000247	C21H26N7O14P2		c	s_1198	-1
% NADH[c]	NADH    bigg.metabolite/nadh;chebi/CHEBI:57945;kegg.compound/C00004;metanetx.chemical/MNXM10;sbo/SBO:0000247	C21H27N7O14P2		c	s_1203	-2
% H+[c]	H+		bigg.metabolite/h;chebi/CHEBI:24636;kegg.compound/C00080;metanetx.chemical/MNXM1;sbo/SBO:0000247	H		c	s_0794	1

%% 3mmol (g CDW)−1s
RD=3;
% RD=0;

% Updated the adjustment of degree of reduciton /Gustav

%Try NADPH for DR balance instead. Also, try to convert NADH to NADPH at an adjustable ratio
% is not required. But we can consider if NADH or NADPH should be used to
% balance the degree of reduction of the biomass
NADH_NADPH=0;

%% NADPH
model.S(find(strcmp(model.mets,'s_1212')),strcmp(model.rxns,'r_4041'))=model.S(find(strcmp(model.mets,'s_1212')),strcmp(model.rxns,'r_4041'))-RD-NADH_NADPH;
%% NADP
 model.S(find(strcmp(model.mets,'s_1207')),strcmp(model.rxns,'r_4041'))=model.S(find(strcmp(model.mets,'s_1207')),strcmp(model.rxns,'r_4041'))+RD+NADH_NADPH;
%% NADH
model.S(find(strcmp(model.mets,'s_1203')),strcmp(model.rxns,'r_4041'))=model.S(find(strcmp(model.mets,'s_1203')),strcmp(model.rxns,'r_4041'))+NADH_NADPH;
%% NAD
 model.S(find(strcmp(model.mets,'s_1198')),strcmp(model.rxns,'r_4041'))=model.S(find(strcmp(model.mets,'s_1198')),strcmp(model.rxns,'r_4041'))-NADH_NADPH;
%% H+
 model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_4041'))=model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_4041'))-RD;


%% Adjust Mw of biomass to 1000

% Biomass_MW=computeMetFormulae(model,'metMwRange','s_0450','fillMets','none','printLevel',0);
% model.S(:,strcmp(model.rxns,'r_4041')) = model.S(:,strcmp(model.rxns,'r_4041'))*1000/mean(Biomass_MW);
% model.S(find(strcmp(model.mets,'s_0450')),strcmp(model.rxns,'r_4041')) = 1;
% 

end

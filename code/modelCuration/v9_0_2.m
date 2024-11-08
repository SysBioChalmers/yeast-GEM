% This scripts applies curations to be applied on yeast-GEM release 9.0.1, to
% get to yeast-GEM release 9.0.2.
% Indicate which Issue/PR are addressed. If multiple curations are performed
% before a new release is made, just add the required code to this script. If
% more extensive coding is required, you can write a separate (generic) function
% that can be kept in the /code/modelCuration folder. Otherwise, try to use
% existing functions whenever possible. In particular /code/curateMetsRxnsGenes
% can do many types of curation.

%% Load yeast-GEM 9.0.1 (requires local yeast-GEM git repository)
model = readYAMLmodel('../../model/yeast-GEM.yml');
% The above comment is temporary, to be used during development. When PR is made, use below code instead.
cd ..
codeDir=pwd();
% model = getEarlierModelVersion('9.0.1');
% model.id='yeastGEM_develop';
% model.version='';
% %dataDir=fullfile(pwd(),'..','data','modelCuration','v9.0.1'); % No dataDir required for these curations
cd modelCuration

%=========================================================================
% We blocked MDH2 in anaerobic conditions (see details in the anerobicModel script)
% Experiments suggest that AKG needs to produced inside the mitochondria
% and exported with ODC1/2 the help or YHM2. After inspecting the kinetics parameters for OAC1, DIC1, YHM2 and
% ODC1/ODC2 (SFC1 is strongly repressed by glucose) we added sulphate and malate as substrate for the OAC1 transporter.
    
%% sulphate[c]	sulphate		bigg.metabolite/so4;chebi/CHEBI:16189;kegg.compound/C00059;metanetx.chemical/MNXM58;sbo/SBO:0000247	O4S		c	s_1467	-2
%  Add sulphate[m]
model = addMetabolite(model, 's_xxxxx', 'sulphate');
model.metComps(end)=9;

%% A proper identifier needs to be added for the metabolite and reaction
% oxaloacetate[c]	oxaloacetate		c	s_1271	-2
% oxaloacetate[m]	oxaloacetate		m	s_1273	-2
% sulfate[m] m	s_xxxxx	-1
% sulfate[c] c	s_1467	-1
model = addReaction(model,'rxxS','metaboliteList',{'s_1271','s_xxxxx',...
    's_1273','s_1467'},'stoichCoeffList',[-1  -1 1 1], 'reversible',true);
model.rev(end)=1;

% (S)-malate[c]	(S)-malate	c	s_0066	-2
% (S)-malate[m]	(S)-malate	m	s_0068	-2
% sulfate[m] m	s_xxxxx	-1
% sulfate[c] c	s_1467	-1
model = addReaction(model,'rxxS2','metaboliteList',{'s_0066','s_xxxxx',...
    's_0068','s_1467'},'stoichCoeffList',[-1  -1 1 1], 'reversible',true);
model.rev(end)=1;


%=========================================================================
% look for all proton symport/antiport reactions and make sure that they
% only enter the cell.

%'s_0794' H+ cytosol
%'s_0796' H+ extracellular
symporterIDs = intersect(find(model.S(find(strcmp(model.mets,'s_0796')),:)),find(model.S(find(strcmp(model.mets,'s_0794')),:)),'stable');
for i = 1:length(symporterIDs)
    if model.S(find(strcmp(model.mets,'s_0796')),symporterIDs(i))<0
        model.lb(symporterIDs(i))=0;
    else
        model.ub(symporterIDs(i))=0;
    end
    if find(strcmp(model.rxns,'r_1258'))==symporterIDs(i) % ignore the sodium transporter, without it, the model does not work
        model.lb(symporterIDs(i))=-1000;
        model.ub(symporterIDs(i))=1000;
    end
end

%=========================================================================
% This section balances reactions and ensures that a correct molecular
% weight can be calculated for the biomass

% Set the charge of all biomass components to 0, should be applied to *all* models, not just anaerobic.
model.metCharges(strcmp(model.mets,'s_3717'))=0; % Protein
model.metCharges(strcmp(model.mets,'s_3718'))=0; % Carbohydrate
model.metCharges(strcmp(model.mets,'s_3719'))=0; % RNA
model.metCharges(strcmp(model.mets,'s_3720'))=0; % DNA
model.metCharges(strcmp(model.mets,'s_3746'))=0; % Lipid backbone
model.metCharges(strcmp(model.mets,'s_3747'))=0; % Lipid chain
model.metCharges(strcmp(model.mets,'s_4205'))=0; % Cofactor
model.metCharges(strcmp(model.mets,'s_4206'))=0; % Ion

% Make the charge of K and Na 1+, should be applied to *all* models, not just anaerobic.
model.metCharges(strcmp(model.mets,'s_1373'))=1;
model.metCharges(strcmp(model.mets,'s_1374'))=1;
model.metCharges(strcmp(model.mets,'s_3776'))=1;
model.metCharges(strcmp(model.mets,'s_1437'))=1;
model.metCharges(strcmp(model.mets,'s_1438'))=1;
model.metCharges(strcmp(model.mets,'s_3775'))=1;

% Balance the charge of all biomass component pseudo reactions by adding the required amount of H+
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_4047')) = -sum(model.S(:,strcmp(model.rxns,'r_4047')).*model.metCharges,'omitnan'); % Protein
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_4049')) = -sum(model.S(:,strcmp(model.rxns,'r_4049')).*model.metCharges,'omitnan'); % RNA
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_4050')) = -sum(model.S(:,strcmp(model.rxns,'r_4050')).*model.metCharges,'omitnan'); % DNA
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_4598')) = -sum(model.S(:,strcmp(model.rxns,'r_4598')).*model.metCharges,'omitnan'); % Cofactor
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_4599')) = -sum(model.S(:,strcmp(model.rxns,'r_4599')).*model.metCharges,'omitnan'); % Ion

% Special case for SLIME rxns
model.metCharges(find(contains(model.metNames,'chain')+contains(model.metNames,'backbone'))) = 0;

% Now, based on the charge balance, find all the reactions that are
% imbalanced, add or remove hydrogen as necessary

% Balance the charge of all imbalanced SLIME reactions by adding the required amount of H+,
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_3975')) = -sum(model.S(:,strcmp(model.rxns,'r_3975')).*model.metCharges,'omitnan'); % Protein
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_3976')) = -sum(model.S(:,strcmp(model.rxns,'r_3976')).*model.metCharges,'omitnan'); % Protein
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_3977')) = -sum(model.S(:,strcmp(model.rxns,'r_3977')).*model.metCharges,'omitnan'); % Protein
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_3978')) = -sum(model.S(:,strcmp(model.rxns,'r_3978')).*model.metCharges,'omitnan'); % Protein
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_4076')) = -sum(model.S(:,strcmp(model.rxns,'r_4076')).*model.metCharges,'omitnan'); % Protein
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_4077')) = -sum(model.S(:,strcmp(model.rxns,'r_4077')).*model.metCharges,'omitnan'); % Protein
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_4078')) = -sum(model.S(:,strcmp(model.rxns,'r_4078')).*model.metCharges,'omitnan'); % Protein
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_4079')) = -sum(model.S(:,strcmp(model.rxns,'r_4079')).*model.metCharges,'omitnan'); % Protein

% Make the charge of HS (hydrogen sulfide) -1
model.metCharges(strcmp(model.mets,'s_0841'))=-1;
model.metCharges(strcmp(model.mets,'s_3906'))=-1;
model.metCharges(strcmp(model.mets,'s_4263'))=-1;

% Balance the reaction r_4702, 'L-cysteine:2-oxoglutarate aminotransferase'
% by adding a proton as reactant
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_4702'))=-1;

% Balance the reaction r_4703, 'L-cysteine:2-oxoglutarate aminotransferase'
% by adding a proton as product
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_4703'))=1;

% Balance the reactions 'r_0774' and 'r_0775', 'NAPRtase' by removing H+ consumption
% and adding a H2O as a reactant
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_0774'))=0;
model.S(find(strcmp(model.mets,'s_0803')),strcmp(model.rxns,'r_0774'))=-1;
model.S(find(strcmp(model.mets,'s_0794')),strcmp(model.rxns,'r_0775'))=0;
model.S(find(strcmp(model.mets,'s_0803')),strcmp(model.rxns,'r_0775'))=-1;

%=========================================================================
% This section focuses on individual reactions that have the wrong
% reversibility/direction/cofactor or should be completley removed

% rename r_0227, it is the plasma membrane ATPase, not a cytosolic ATPase
model.rxnNames(strcmp(model.rxns,'r_0227')) = {'ATPase, plasma membrane'};

% make sure both formate-THF ligases are reversible.
model.lb(strcmp(model.rxns,'r_0446')) = -1000;

% ADE3 and MIS1, methylenetetrahydrofolate dehydrogenase (NADP+)  [EC 1.5.1.5]
% make irrevrersible
model.ub(strcmp(model.rxns,'r_0732')) = 0;
model.ub(strcmp(model.rxns,'r_0733')) = 0;

% There is no evidence for this PFK1 side reaction in yeast. Consider
% removing the reaction completley
model.ub(strcmp(model.rxns,'r_0887')) = 0;

% TYR1 incorrectly annotated as using NAD, should be NADP
model.S(find(strcmp(model.mets,'s_1212')),strcmp(model.rxns,'r_0939'))=0; %NADPH
model.S(find(strcmp(model.mets,'s_1207')),strcmp(model.rxns,'r_0939'))=0; %NADP
model.S(find(strcmp(model.mets,'s_1203')),strcmp(model.rxns,'r_0939'))=1; %NADH
model.S(find(strcmp(model.mets,'s_1198')),strcmp(model.rxns,'r_0939'))=-1;%NAD

% Make esterification reactions irreversible. positive deltaG
model.ub(strcmp(model.rxns,'r_4713')) = 0; %diethyl succinate
model.ub(strcmp(model.rxns,'r_4714')) = 0; %monoethyl succinate

% Make polyphosphate hydrolase and diphosphate transport over cell membrane
% both irreversible
model = setParam(model,'lb',{'r_4723','r_4724','r_4725'},0);
model = setParam(model,'lb','r_4460',0);

% In principle only ammonium (not ammonia) can be transported tomitochondria. In the cytosol, NH4+ is the main form. 
% Additionally, NH3 transport processes between different compartments are assumed to operate close to thermodynamic equilibrium
% -and since no transport proteins that could translocate NHX between compartments are described in literature, passive diffusion 
% of NH3 between vacuole and cytosol, as well as between cytosol and mitochondria, are assumed.
% Cueto-Rojas, Hugo F., et al. "Membrane potential independent transport of NH 3 in the absence of ammonium permeases in Saccharomyces cerevisiae." BMC Systems Biology 11 (2017): 1-13. 
%'r_1965'	'NH3 transport'	'ammonium[cytoplasm]  <=> ammonium[mitochondrion] '	0.012058369	-0.277342481	'' 
model = setParam(model,'eq',{'r_1965'},0);


%=========================================================================
%% Condition-specific gene expression. These can be enabled with scripts

% % Glycine cleavage only active when glycine is used as nitrogen source
% model.ub(strcmp(model.rxns,'r_0501'))=0; %glycine cleavage, mitochondrion
% model.lb(strcmp(model.rxns,'r_0501'))=0;
% model.ub(strcmp(model.rxns,'r_0507'))=0; %glycine cleavage complex (lipoylprotein), mitochondrion
% model.lb(strcmp(model.rxns,'r_0507'))=0;
% model.ub(strcmp(model.rxns,'r_0509'))=0; %glycine cleavage complex (lipoamide), mitochondrion
% model.lb(strcmp(model.rxns,'r_0509'))=0;

% Glutamate synthase repressed in excess nitrogen
model.ub(strcmp(model.rxns,'r_0472'))=0;

% The carnitine shuttle requires exogeneous carnitine
% r_0252	'carnitine O-acetyltransferase'	0	1000	'(R)-carnitine[cytoplasm] + acetyl-CoA[cytoplasm]  -> coenzyme A[cytoplasm] + O-acetylcarnitine[cytoplasm] '	0	0	0	0	0.167350072
model = setParam(model,'eq',{'r_0252'},0);

%=========================================================================

%% DO NOT CHANGE OR REMOVE THE CODE BELOW THIS LINE.
% Show some metrics:
% cd(fullfile(codeDir,'modelTests'))
% disp('Run gene essentiality analysis')
% [new.accuracy,new.tp,new.tn,new.fn,new.fp] = essentialGenes(model);
% fprintf('Genes in model: %d\n',numel(model.genes));
% fprintf('Gene essentiality accuracy: %.4f\n', new.accuracy);
% fprintf('True non-essential genes: %d\n', numel(new.tp));
% fprintf('True essential genes: %d\n', numel(new.tn));
% fprintf('False non-essential genes: %d\n', numel(new.fp));
% fprintf('False essential genes: %d\n', numel(new.fn));
% fprintf('\nRun growth analysis\n')
% R2=growth(model);
% fprintf('R2 of growth prediction: %.4f\n', R2);
%
% Save model:
% cd ..
% saveYeastModel(model)
% cd modelCuration

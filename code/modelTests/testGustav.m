clear;
close all;
%%
model=loadYeastModel;
% cd ../modelCuration/
% v9_1_0;
% cd ../modelTests/
% Convert to anaerobic
% model = anaerobicModel(model);

model = anaerobicModel_new(model,false);

% Glucose uptake rate
model = setParam(model,'eq','r_1714',-23);

% Solve
res=solveLP(model,1);
FLUX=res.x;

v_glc=FLUX(getIndexes(model,'r_1714','rxns'),:);
v_eth=FLUX(getIndexes(model,'r_1761','rxns'),:);
v_CO2=FLUX(getIndexes(model,'r_1672','rxns'),:);
v_gly=FLUX(getIndexes(model,'r_1808','rxns'),:);
v_growth=FLUX(getIndexes(model,'r_4041','rxns'),:);
%%
%glycerol ethanol Co2
%4.5 ± 0.4  31 ± 2  38 ± 10
data=[4.5 31*1.2 38 0.36];
sim=[v_gly v_eth v_CO2 v_growth];
errorVal=[0.4 2 10 0.02];
b1=bar(data./data,'FaceAlpha',0.5);hold on;b2=bar(sim./data,'FaceAlpha',0.5);

hold on

er = errorbar([1 2 3 4],data./data,errorVal./data,errorVal./data);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  


legend({'data','simulation'});
ylabel('Relative value');
xticklabels({'Glycerol','Ethanol','CO2','Biomass'})
% 
% temp_model=model;
% temp_model.metNames=strcat( model.metNames, repmat('[',length(model.mets),1),model.compNames(temp_model.metComps),repmat(']',length(model.mets),1) ); 
% 
% 
% % Pack everything into a nice table
% rxns_reacs=printRxnFormula(temp_model,'rxnAbbrList',model.rxns,'metNameFlag',1,'printFlag',0);
% [massImbalance, imBalancedMass, imBalancedCharge, imBalancedRxnBool, elements, missingFormulaeBool, balancedMetBool] = checkMassChargeBalance(model);
% tab=table(model.rxns,model.rxnNames,rxns_reacs,abs(FLUX./v_glc),FLUX,model.grRules,imBalancedRxnBool,imBalancedCharge,imBalancedMass);
% %tab=table(model.rxns,model.rxnNames,rxns_reacs,abs(FLUX./v_glc),FLUX,model.grRules);
% tab = sortrows(tab,"Var4","descend");
% 
% printFluxes(model,FLUX,true)
% 
% 
% 
% 
% 
% %% Calculate formula and degree of reduciton of biomass
% 
[mwRange,metFormulae,elements,metEle]=computeMetFormulae(model,'metMwRange','s_0450','fillMets','none','printLevel',0);
mwRange

% Print out the results of the 
Biomass_index = find(strcmp(model.metNames,'biomass'));

%Degree of reduction per element. Order of the elements 'C', 'H', 'O', 'N'
DR_per_ele = [4, 1, -2, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
DR_per_Cmol=sum(metEle(Biomass_index,:).*DR_per_ele)/metEle(Biomass_index,1)

Biomass_formula_Cmol=metEle(Biomass_index,:)/metEle(Biomass_index,1);

MW_per_Cmol_min = mwRange/metEle(Biomass_index,1)


function model_an=anaerobicModel_new (model,NLIM)
    model_an=model;
    NLIM=false;
    modification_level=4;
    
    if modification_level>0
    
        % this script makes corrections that should be applied to *all* models, not
        % just the anaerobic. Adjustments to the anaerobic model are found in the
        % anaerobicModel function
    
    
        %=========================================================================
    
        % Create a pathway for anaplerosis. Check alternatives? this one might work
        % on other diacids also? If considered important for "flux reasons", dig
        % into this transporter literature
        % 'r_1264' 'succinate transport'   'phosphate[mitochondrion] + succinate[cytoplasm]  -> phosphate[cytoplasm] + succinate[mitochondrion] '
        model.lb(findRxnIDs(model,'r_1264'))=-1000; %DIC1. I am not sure its actually a phosphate/succinate antiport. It not clear from the article. It could also be specific for both phosphate and succinate
    
        %=========================================================================
    
        % 
        % % look for all proton symport/antiport reactions and make sure that they
        % % only enter the cell.
        % 
        % %'s_0794' H+ cytosol
        % %'s_0796' H+ extracellular
        symporterIDs = intersect(find(model_an.S(find(strcmp(model_an.mets,'s_0796')),:)),find(model_an.S(find(strcmp(model_an.mets,'s_0794')),:)),'stable');
        for i = 1:length(symporterIDs)
            if model_an.S(find(strcmp(model_an.mets,'s_0796')),symporterIDs(i))<0
                model_an.lb(symporterIDs(i))=0;
            else
                model_an.ub(symporterIDs(i))=0;
            end

            if find(strcmp(model_an.rxns,'r_1258'))==symporterIDs(i) % ignore the sodium transporter, without it, the model doesn't work
                model_an.lb(symporterIDs(i))=-1000;
                model_an.ub(symporterIDs(i))=1000;
            end
        end
        % 
        %=========================================================================
    
        % This section balances reactions and ensures that a correct molecular
        % weight can be calculated for the biomass
        % 
        % Set the charge of all biomass components to 0, should be applied to *all* models, not just anaerobic.
        model_an.metCharges(strcmp(model_an.mets,'s_3717'))=0; % Protein
        model_an.metCharges(strcmp(model_an.mets,'s_3718'))=0; % Carbohydrate
        model_an.metCharges(strcmp(model_an.mets,'s_3719'))=0; % RNA
        model_an.metCharges(strcmp(model_an.mets,'s_3720'))=0; % DNA
        model_an.metCharges(strcmp(model_an.mets,'s_3746'))=0; % Lipid backbone
        model_an.metCharges(strcmp(model_an.mets,'s_3747'))=0; % Lipid chain
        model_an.metCharges(strcmp(model_an.mets,'s_4205'))=0; % Cofactor
        model_an.metCharges(strcmp(model_an.mets,'s_4206'))=0; % Ion

        % Make the charge of K and Na 1+, should be applied to *all* models, not just anaerobic.
        model_an.metCharges(strcmp(model_an.mets,'s_1373'))=1;
        model_an.metCharges(strcmp(model_an.mets,'s_1374'))=1;
        model_an.metCharges(strcmp(model_an.mets,'s_3776'))=1;
        model_an.metCharges(strcmp(model_an.mets,'s_1437'))=1;
        model_an.metCharges(strcmp(model_an.mets,'s_1438'))=1;
        model_an.metCharges(strcmp(model_an.mets,'s_3775'))=1;

        % Balance the charge of all biomass component pseudo reactions by adding the required amount of H+
        model_an.S(find(strcmp(model_an.mets,'s_0794')),strcmp(model_an.rxns,'r_4047')) = -sum(model_an.S(:,strcmp(model_an.rxns,'r_4047')).*model_an.metCharges,'omitnan'); % Protein
        model_an.S(find(strcmp(model_an.mets,'s_0794')),strcmp(model_an.rxns,'r_4049')) = -sum(model_an.S(:,strcmp(model_an.rxns,'r_4049')).*model_an.metCharges,'omitnan'); % RNA
        model_an.S(find(strcmp(model_an.mets,'s_0794')),strcmp(model_an.rxns,'r_4050')) = -sum(model_an.S(:,strcmp(model_an.rxns,'r_4050')).*model_an.metCharges,'omitnan'); % DNA
        model_an.S(find(strcmp(model_an.mets,'s_0794')),strcmp(model_an.rxns,'r_4598')) = -sum(model_an.S(:,strcmp(model_an.rxns,'r_4598')).*model_an.metCharges,'omitnan'); % Cofactor
        model_an.S(find(strcmp(model_an.mets,'s_0794')),strcmp(model_an.rxns,'r_4599')) = -sum(model_an.S(:,strcmp(model_an.rxns,'r_4599')).*model_an.metCharges,'omitnan'); % Ion


        % Special case for SLIME rxns
        model_an.metCharges(find(contains(model_an.metNames,'chain')+contains(model_an.metNames,'backbone'))) = 0;
    
        % Now, based on the charge balance, find all the reactions that are
        % imbalanced, add or remove hydrogen as necessary
    

        % Balance the charge of all imbalanced SLIME reactions by adding the required amount of H+,
        model_an.S(find(strcmp(model_an.mets,'s_0794')),strcmp(model_an.rxns,'r_3975')) = -sum(model_an.S(:,strcmp(model_an.rxns,'r_3975')).*model_an.metCharges,'omitnan'); % Protein
        model_an.S(find(strcmp(model_an.mets,'s_0794')),strcmp(model_an.rxns,'r_3976')) = -sum(model_an.S(:,strcmp(model_an.rxns,'r_3976')).*model_an.metCharges,'omitnan'); % Protein
        model_an.S(find(strcmp(model_an.mets,'s_0794')),strcmp(model_an.rxns,'r_3977')) = -sum(model_an.S(:,strcmp(model_an.rxns,'r_3977')).*model_an.metCharges,'omitnan'); % Protein
        model_an.S(find(strcmp(model_an.mets,'s_0794')),strcmp(model_an.rxns,'r_3978')) = -sum(model_an.S(:,strcmp(model_an.rxns,'r_3978')).*model_an.metCharges,'omitnan'); % Protein
        model_an.S(find(strcmp(model_an.mets,'s_0794')),strcmp(model_an.rxns,'r_4076')) = -sum(model_an.S(:,strcmp(model_an.rxns,'r_4076')).*model_an.metCharges,'omitnan'); % Protein
        model_an.S(find(strcmp(model_an.mets,'s_0794')),strcmp(model_an.rxns,'r_4077')) = -sum(model_an.S(:,strcmp(model_an.rxns,'r_4077')).*model_an.metCharges,'omitnan'); % Protein
        model_an.S(find(strcmp(model_an.mets,'s_0794')),strcmp(model_an.rxns,'r_4078')) = -sum(model_an.S(:,strcmp(model_an.rxns,'r_4078')).*model_an.metCharges,'omitnan'); % Protein
        model_an.S(find(strcmp(model_an.mets,'s_0794')),strcmp(model_an.rxns,'r_4079')) = -sum(model_an.S(:,strcmp(model_an.rxns,'r_4079')).*model_an.metCharges,'omitnan'); % Protein


        % Make the charge of HS (hydrogen sulfide) -1
        model_an.metCharges(strcmp(model_an.mets,'s_0841'))=-1;
        model_an.metCharges(strcmp(model_an.mets,'s_3906'))=-1;
        model_an.metCharges(strcmp(model_an.mets,'s_4263'))=-1;

        % Balance the reaction r_4702, 'L-cysteine:2-oxoglutarate aminotransferase'
        % by adding a proton as reactant
        model_an.S(find(strcmp(model_an.mets,'s_0794')),strcmp(model_an.rxns,'r_4702'))=-1;

        % Balance the reaction r_4703, 'L-cysteine:2-oxoglutarate aminotransferase'
        % by adding a proton as product
        model_an.S(find(strcmp(model_an.mets,'s_0794')),strcmp(model_an.rxns,'r_4703'))=1;

        % Balance the reactions 'r_0774' and 'r_0775', 'NAPRtase' by removing H+ consumption
        % and adding a H2O as a reactant
        model_an.S(find(strcmp(model_an.mets,'s_0794')),strcmp(model_an.rxns,'r_0774'))=0;
        model_an.S(find(strcmp(model_an.mets,'s_0803')),strcmp(model_an.rxns,'r_0774'))=-1;
        model_an.S(find(strcmp(model_an.mets,'s_0794')),strcmp(model_an.rxns,'r_0775'))=0;
        model_an.S(find(strcmp(model_an.mets,'s_0803')),strcmp(model_an.rxns,'r_0775'))=-1;
    end
    
    if modification_level>1
        %=========================================================================
        % % This section focuses on individual reactions that have the wrong
        % % reversibility/direction/cofactor or should be completley removed
        % 
        % % rename r_0227, it is the plasma membrane ATPase, not a cytosolic ATPase
        % model_an.rxnNames(strcmp(model_an.rxns,'r_0227')) = {'ATPase, plasma membrane'};
        % 
        % % make sure both formate-THF ligases are reversible.
        % model_an.lb(strcmp(model_an.rxns,'r_0446')) = -1000;
        % 
        % % ADE3 and MIS1, methylenetetrahydrofolate dehydrogenase (NADP+)  [EC 1.5.1.5]
        % % make irrevrersible
        % model_an.ub(strcmp(model_an.rxns,'r_0732')) = 0;
        % model_an.ub(strcmp(model_an.rxns,'r_0733')) = 0;
        % 
        % % There is no evidence for this PFK1 side reaction in yeast. Consider
        % % removing the reaction completley
        % model_an.ub(strcmp(model_an.rxns,'r_0887')) = 0;
        % 
        % % TYR1 incorrectly annotated as using NAD, should be NADP
        % model_an.S(find(strcmp(model_an.mets,'s_1212')),strcmp(model_an.rxns,'r_0939'))=0; %NADPH
        % model_an.S(find(strcmp(model_an.mets,'s_1207')),strcmp(model_an.rxns,'r_0939'))=0; %NADP
        % model_an.S(find(strcmp(model_an.mets,'s_1203')),strcmp(model_an.rxns,'r_0939'))=1; %NADH
        % model_an.S(find(strcmp(model_an.mets,'s_1198')),strcmp(model_an.rxns,'r_0939'))=-1;%NAD
        % 
        % % Make esterification reactions irreversible. positive deltaG
        % model_an.ub(strcmp(model_an.rxns,'r_4713')) = 0; %diethyl succinate
        % model_an.ub(strcmp(model_an.rxns,'r_4714')) = 0; %monoethyl succinate
        % 
        % % Make polyphosphate hydrolase and diphosphate transport over cell membrane
        % % both irreversible
        % model_an = setParam(model_an,'lb',{'r_4723','r_4724','r_4725'},0);
        % model_an = setParam(model_an,'lb','r_4460',0);
        % %=========================================================================
    end
    
    % 
    % GAM   = 55;
    % P     = 0.461;  %Data from Nissen et al. 1997
    % NGAM  = 1;
    % 
    % model_an = changeGAM(model_an,GAM,NGAM);
    % if(NLIM)
    %     model_an = scaleBioMass(model_an,'protein',0.289,'',false);
    %     model_an = scaleBioMass(model_an,'lipid',0.048,'',false);
    %     model_an = scaleBioMass(model_an,'RNA',0.077,'carbohydrate',false);
    % 
    % else        
    %     model_an = scaleBioMass(model_an,'protein',P,'carbohydrate',false);
    % 
    % end
    
    % %2nd change: Removes the requirement of heme a, NAD(PH), coenzyme A in the biomass equation
    % %            (not used under anaerobic conditions)
    
    % I disagree. Please do not remove NAD(P)H and CoA, include the exhange
    % reactions for the relevant vitamins instead. See below /Gustav
    mets = {'s_3714'}; %,'s_1198','s_1203','s_1207','s_1212','s_0529'};
    
    [~,met_index] = ismember(mets,model_an.mets);
    model_an.S(met_index,strcmp(model_an.rxns,'r_4598')) = 0;
    
    %3rd change: Changes media to anaerobic (no O2 uptake and allows sterol
    %            and fatty acid exchanges)
    disp('bla')
    model_an.lb(strcmp(model_an.rxns,'r_1992')) = 0;        %O2
    model_an.lb(strcmp(model_an.rxns,'r_1757')) = -1000;    %ergosterol
    model_an.lb(strcmp(model_an.rxns,'r_1915')) = -1000;    %lanosterol
    model_an.lb(strcmp(model_an.rxns,'r_1994')) = -1000;    %palmitoleate
    model_an.lb(strcmp(model_an.rxns,'r_2106')) = -1000;    %zymosterol
    model_an.lb(strcmp(model_an.rxns,'r_2134')) = -1000;    %14-demethyllanosterol
    %remove this due to NADH recycling to ergosterol
    model_an.lb(strcmp(model_an.rxns,'r_2137')) = 0;    %ergosta-5,7,22,24(28)-tetraen-3beta-ol
    model_an.lb(strcmp(model_an.rxns,'r_2189')) = -1000;    %oleate
    
    % Added exchange of vitamins enabling NAD(P)H and CoA syntheis in anaerobic
    % conditions /Gustav
    model_an.lb(strcmp(model_an.rxns,'r_1967')) = -1000;    %nicotinate
    model_an.lb(strcmp(model_an.rxns,'r_1548')) = -1000;    %(R)-pantothenate
    
    
    if modification_level>2
        
        %% Changes to the model that give correct phenotype for anaerobic batch growht on minimal glucose media
        % Inhibit MDH2 during excess glucose = anaerobic conditions.
        model_an = setParam(model_an,'lb','r_0714',0);
        model_an = setParam(model_an,'ub','r_0714',0);
    
    
    
        % GCY1 has a positive DeltaG and is part of a transhydrogenase cycle NADH -> NADPH
        model_an.ub(strcmp(model_an.rxns,'r_0487')) = 0;
        
        %Glutamate synthase repressed in excess nitrogen
        model_an.ub(strcmp(model_an.rxns,'r_0472'))=0;
    
        %glycine cleavage respressed by presence of glucose
        model_an.ub(strcmp(model_an.rxns,'r_0501'))=0; %glycine cleavage, mitochondrion
        model_an.lb(strcmp(model_an.rxns,'r_0501'))=0;
        model_an.ub(strcmp(model_an.rxns,'r_0507'))=0; %glycine cleavage complex (lipoylprotein), mitochondrion
        model_an.lb(strcmp(model_an.rxns,'r_0507'))=0;
        model_an.ub(strcmp(model_an.rxns,'r_0509'))=0; %glycine cleavage complex (lipoamide), mitochondrion
        model_an.lb(strcmp(model_an.rxns,'r_0509'))=0;
    
        % MAE1 and IDP are likely not major mitochondrial NADPH sources.
        %model_an.ub(strcmp(model_an.rxns,'r_0719')) =0; % malic enzyme (MAE1), mitochondrion
        %model_an.ub(strcmp(model_an.rxns,'r_2131')) = 0; % isocitrate dehydrogenase (IDP1), mitochondrion
    
    end
    
    if modification_level>3
        model_an.lb(strcmp(model_an.rxns,'r_0226')) = -1000; % ATP synthase
    
        
        %% 'r_0659'	'isocitrate dehydrogenase (NADP)'	'isocitrate[cytoplasm] + NADP(+)[cytoplasm]  <=> 2-oxoglutarate[cytoplasm] + carbon dioxide[cytoplasm] + NADPH[cytoplasm] '	0.014152441	0.325506132	'YLR174W'
        model_an = setParam(model_an,'eq',{'r_0659'},0);
    
        %% r_0252	'carnitine O-acetyltransferase'	0	1000	'(R)-carnitine[cytoplasm] + acetyl-CoA[cytoplasm]  -> coenzyme A[cytoplasm] + O-acetylcarnitine[cytoplasm] '	0	0	0	0	0.167350072
        model_an = setParam(model_an,'eq',{'r_0252'},0);
        
        %% sulphate[c]	sulphate		bigg.metabolite/so4;chebi/CHEBI:16189;kegg.compound/C00059;metanetx.chemical/MNXM58;sbo/SBO:0000247	O4S		c	s_1467	-2
        model_an = addMetabolite(model_an, 's_xxxxx', 'sulphate');
        model_an.metComps(end)=9;
        
        % oxaloacetate[c]	oxaloacetate		c	s_1271	-2
        % oxaloacetate[m]	oxaloacetate		m	s_1273	-2
        % sulfate[m] m	s_xxxxx	-1
        % sulfate[c] c	s_1467	-1
        model_an = addReaction(model_an,'rxxS','metaboliteList',{'s_1271','s_xxxxx',...
            's_1273','s_1467'},'stoichCoeffList',[-1  -1 1 1], 'reversible',true);
        model_an.rev(end)=1;
    
        % (S)-malate[c]	(S)-malate	c	s_0066	-2
        % (S)-malate[m]	(S)-malate	m	s_0068	-2
        % sulfate[m] m	s_xxxxx	-1
        % sulfate[c] c	s_1467	-1
        model_an = addReaction(model_an,'rxxS2','metaboliteList',{'s_0066','s_xxxxx',...
            's_0068','s_1467'},'stoichCoeffList',[-1  -1 1 1], 'reversible',true);
        model_an.rev(end)=1;
    
    
    end
    
    %% A reaction converting NADH to NAD + at 3 mmol (g CDW)−1 was coupled to the growth reaction to give the correct ratio of glycerol production to
    % glucose consumption, and to align the degree of reduction (as defined by (Heijnen, 1994)) of the modeled biomass to a published value for
    % degree of reduction of S. cerevisiae biomass at 4.2 C-mol−1 (Lange and Heijnen, 2001).
    % NAD[c]	NAD		bigg.metabolite/nad;chebi/CHEBI:57540;kegg.compound/C00003;metanetx.chemical/MNXM8;sbo/SBO:0000247	C21H26N7O14P2		c	s_1198	-1
    % NADH[c]	NADH    bigg.metabolite/nadh;chebi/CHEBI:57945;kegg.compound/C00004;metanetx.chemical/MNXM10;sbo/SBO:0000247	C21H27N7O14P2		c	s_1203	-2
    % H+[c]	H+		bigg.metabolite/h;chebi/CHEBI:24636;kegg.compound/C00080;metanetx.chemical/MNXM1;sbo/SBO:0000247	H		c	s_0794	1
    
    % % 3mmol (g CDW)−1s
    % RD=3;
    % 
    % % Updated the adjustment of degree of reduciton /Gustav
    % % Try NADPH for DR balance instead. Also, try to convert NADH to NADPH at an adjustable ratio
    % % is not required. But we can consider if NADH or NADPH should be used to
    % % balance the degree of reduction of the biomass
    % NADH_NADPH=0;
    % 
    % % NADPH
    % model_an.S(find(strcmp(model_an.mets,'s_1212')),strcmp(model_an.rxns,'r_4041'))=model_an.S(find(strcmp(model_an.mets,'s_1212')),strcmp(model_an.rxns,'r_4041'))-RD-NADH_NADPH;
    % % NADP
    % model_an.S(find(strcmp(model_an.mets,'s_1207')),strcmp(model_an.rxns,'r_4041'))=model_an.S(find(strcmp(model_an.mets,'s_1207')),strcmp(model_an.rxns,'r_4041'))+RD+NADH_NADPH;
    % % NADH
    % model_an.S(find(strcmp(model_an.mets,'s_1203')),strcmp(model_an.rxns,'r_4041'))=model_an.S(find(strcmp(model_an.mets,'s_1203')),strcmp(model_an.rxns,'r_4041'))+NADH_NADPH;
    % % NAD
    % model_an.S(find(strcmp(model_an.mets,'s_1198')),strcmp(model_an.rxns,'r_4041'))=model_an.S(find(strcmp(model_an.mets,'s_1198')),strcmp(model_an.rxns,'r_4041'))-NADH_NADPH;
    % % H+
    % model_an.S(find(strcmp(model_an.mets,'s_0794')),strcmp(model_an.rxns,'r_4041'))=model_an.S(find(strcmp(model_an.mets,'s_0794')),strcmp(model_an.rxns,'r_4041'))-RD;

    %% Adjust Mw of biomass to 1000
    % Biomass_MW=computeMetFormulae(model_an,'metMwRange','s_0450','fillMets','none','printLevel',0);
    % model_an.S(:,strcmp(model_an.rxns,'r_4041')) = model_an.S(:,strcmp(model_an.rxns,'r_4041'))*1000/mean(Biomass_MW);
    % model_an.S(find(strcmp(model_an.mets,'s_0450')),strcmp(model_an.rxns,'r_4041')) = 1;
    % 
    % model_an.lb(strcmp(model_an.rxns,'r_1714')) = -23; % glucose uptake
    % sol=solveLP(model_an);

end
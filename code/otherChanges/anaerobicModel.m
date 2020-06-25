%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% anaerobicModel.m
% Converts model to anaerobic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = anaerobicModel(model)

%1st change: Refit GAM and NGAM to exp. data, change biomass composition
GAM   = 30.49;  %Data from Nissen et al. 1997
P     = 0.461;  %Data from Nissen et al. 1997
NGAM  = 0;      %Refit done in Jouthen et al. 2012

model = changeGAM(model,GAM,NGAM);
model = scaleBioMass(model,'protein',P,'carbohydrate');

%2nd change: Removes the requirement of heme a, NAD(PH), coenzyme A in the biomass equation
%            (not used under anaerobic conditions)
mets = {'s_3714[c]','s_1198[c]','s_1203[c]','s_1207[c]','s_1212[c]','s_0529[c]'};
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
model.lb(strcmp(model.rxns,'r_2137')) = -1000;    %ergosta-5,7,22,24(28)-tetraen-3beta-ol
model.lb(strcmp(model.rxns,'r_2189')) = -1000;    %oleate

%4th change: Blocked pathways for proper glycerol production
%Block oxaloacetate-malate shuttle (not present in anaerobic conditions)
model.lb(strcmp(model.rxns,'r_0713')) = 0; %Mithocondria
model.lb(strcmp(model.rxns,'r_0714')) = 0; %Cytoplasm
%Block glycerol dehydroginase (only acts in microaerobic conditions)
model.ub(strcmp(model.rxns,'r_0487')) = 0;
%Block 2-oxoglutarate + L-glutamine -> 2 L-glutamate (alternative pathway)
model.ub(strcmp(model.rxns,'r_0472')) = 0;

end

%%

function model = changeGAM(model,GAM,NGAM)

bioPos = strcmp(model.rxnNames,'biomass pseudoreaction');
for i = 1:length(model.mets)
    S_ix  = model.S(i,bioPos);
    isGAM = sum(strcmp({'ATP [cytoplasm]','ADP [cytoplasm]','H2O [cytoplasm]', ...
        'H+ [cytoplasm]','phosphate [cytoplasm]'},model.metNames{i})) == 1;
    if S_ix ~= 0 && isGAM
        model.S(i,bioPos) = sign(S_ix)*GAM;
    end
end

if nargin >1
    pos = strcmp(model.rxnNames,'non-growth associated maintenance reaction');%NGAM
    model = setParam(model,'eq',model.rxns(pos),NGAM);% set both lb and ub to be mu
end

end

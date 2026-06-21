function model = anaerobicModelOld(model)
% anaerobicModelOld
%   This file has been replaced with anaerobicModel since yeast-GEM 9.1.0.
%   This function is kept just for comparison purposes.
%
%   Inputs: model           (struct) aerobic model
%   Output: model           (struct) anaerobic model
%   
%   Usage: model = anaerobicModel(model)
%

%1st change: Refit GAM and NGAM to exp. data, change biomass composition
GAM   = 30.49;  %Data from Nissen et al. 1997
P     = 0.461;  %Data from Nissen et al. 1997
NGAM  = 0;      %Refit done in Jouthen et al. 2012

model = changeGAM(model,GAM,NGAM);
model = scaleBioMassOld(model,'protein',P,'carbohydrate',false);

%2nd change: Removes the requirement of heme a, NAD(PH), coenzyme A in the biomass equation
%            (not used under anaerobic conditions)
mets = {'s_3714','s_1198','s_1203','s_1207','s_1212','s_0529'};
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

function model = scaleBioMassOld(model,component,new_value,balance_out,dispOutput)
  % scaleBioMass
  %   Scales the biomass composition
  %
  %   model          (struct) metabolic model in COBRA format
  %   component      (str) name of the component to rescale (e.g. "protein")
  %   new_value      (float) new total fraction for said component
  %   balance_out    (str, opt) if chosen, the name of another component with which
  %                  the model will be balanced out so that the total mass remains = 1 g/gDW
  %                  provide empty string '' if this should not be done
  %   dispOutput     (bool, opt) if output from sumBioMass should be displayed (default = true)
  %
  %   model          (struct) modified model
  %
  %   Usage: model = scaleBioMass(model,component,new_value,balance_out,dispOutput)
  %

if nargin < 5
    dispOutput = true;
end
if nargin < 4
    balance_out = '';
end
  
%Measure current composition and rescale:
[~,P,C,R,D,L,I,F] = sumBioMassOld(model,dispOutput);
content_all = {'carbohydrate','protein','lipid','RNA','DNA','ion','cofactor'};
content_Cap = {'C','P','L','R','D','I','F'};
pos         = strcmp(content_all,component);
old_value   = eval(content_Cap{pos});
f           = new_value / old_value;
model       = rescalePseudoReaction(model,component,f);

%Balance out (if desired):
if ~isempty(balance_out)
    pos           = strcmp(content_all,balance_out);
    balance_value = eval(content_Cap{pos});
    f             = (balance_value - (new_value - old_value)) / balance_value;
    model         = rescalePseudoReaction(model,balance_out,f);
end
end

function [X,P,C,R,D,L,I,F] = sumBioMassOld(model,dispOutput)
  % sumBioMass
  %   Calculates breakdown of biomass
  %
  %   model         (struct) Metabolic model in COBRA format
  %   dispOutput    (bool, opt) If output should be displayed (default = true)
  %
  %   X             (float) Total biomass fraction [gDW/gDW]
  %   P             (float) Protein fraction [g/gDW]
  %   C             (float) Carbohydrate fraction [g/gDW]
  %   R             (float) RNA fraction [g/gDW]
  %   D             (float) DNA fraction [g/gDW]
  %   L             (float) Lipid fraction [g/gDW]
  %   F             (float) cofactor [g/gDW]
  %   I             (float) ion [g/gDW]
  %
  %   Usage: [X,P,C,R,D,L,I,F] = sumBioMass(model,dispOutput)
  %
  %   Function adapted from SLIMEr: https://github.com/SysBioChalmers/SLIMEr
  %

if nargin < 2
    dispOutput = true;
end

%Load original biomass component MWs:
%TODO: compute MW automatically from chemical formulas (check that all components have them first)
fid = fopen('../../data/physiology/biomassComposition_Forster2003.tsv');
Forster2003 = textscan(fid,'%s %s %f32 %f32 %s','Delimiter','\t','HeaderLines',1);
data.mets   = Forster2003{1};
data.MWs    = double(Forster2003{4});
fclose(fid);

%load additional cofactor/ion MWs:
fid = fopen('../../data/physiology/biomassComposition_Cofactor_Ion.tsv');
CofactorsIons = textscan(fid,'%s %s %f32 %f32 %s %s','Delimiter','\t','HeaderLines',1);
data_new.mets = CofactorsIons{1};
data_new.MWs  = double(CofactorsIons{4});
fclose(fid);
for i = 1:length(data_new.mets)
    if ~ismember(data_new.mets(i),data.mets)
        data.mets = [data.mets; data_new.mets(i)];
        data.MWs  = [data.MWs; data_new.MWs(i)];
    end
end

%Get main fractions:
[P,X] = getFraction(model,data,'P',0,dispOutput);
[C,X] = getFraction(model,data,'C',X,dispOutput);
[R,X] = getFraction(model,data,'R',X,dispOutput);
[D,X] = getFraction(model,data,'D',X,dispOutput);
[L,X] = getFraction(model,data,'L',X,dispOutput);
[I,X] = getFraction(model,data,'I',X,dispOutput);
[F,X] = getFraction(model,data,'F',X,dispOutput);

if dispOutput
    disp(['X -> ' num2str(X) ' gDW/gDW'])
    % Simulate growth:
    sol = solveLP(model,1);
    disp(['Growth = ' num2str(sol.f) ' 1/h'])
    disp(' ')
end

end

%%

function [F,X] = getFraction(model,data,compType,X,dispOutput)

%Define pseudoreaction name:
rxnName = [compType ' pseudoreaction'];
rxnName = strrep(rxnName,'P','protein');
rxnName = strrep(rxnName,'C','carbohydrate');
rxnName = strrep(rxnName,'N','biomass');
rxnName = strrep(rxnName,'L','lipid backbone');
rxnName = strrep(rxnName,'R','RNA');
rxnName = strrep(rxnName,'D','DNA');
rxnName = strrep(rxnName,'I','ion');
rxnName = strrep(rxnName,'F','cofactor');

%Add up fraction:
rxnPos = strcmp(model.rxnNames,rxnName);
if ~all(rxnPos==0)
    isSub   = model.S(:,rxnPos) < 0;        %substrates in pseudo-rxn
    if strcmp(compType,'L')
        F = -sum(model.S(isSub,rxnPos));   %g/gDW
    else
        F = 0;
        %Add up all components:
        for i = 1:length(model.mets)
            pos = strcmp(data.mets,model.mets{i});
            if isSub(i) && sum(pos) == 1
                if strcmp(compType,'I') || strcmp(compType,'F')
                    MW = data.MWs(pos);
                else
                    MW = data.MWs(pos)-18;
                end
                abundance = -model.S(i,rxnPos)*MW/1000;
                F         = F + abundance;
            end
        end
    end
    X = X + F;
    
    if dispOutput
        disp([compType ' -> ' num2str(F) ' g/gDW'])
    end
else
    if dispOutput
        disp([compType ' does not exist '])
    end
    F = 0;
    X = X + F;
end
end
clear; close all
% Load the old and new model
model902 = getEarlierModelVersion('9.0.2');
model910 = getEarlierModelVersion('9.1.0');

%% Run growth tests
funcs   = {@growth, @growthOld, @growth, @growthOld};
models  = {model910, model910, model902, model902};
titles  = {'v9.1.0 model, new anaerobic script', ...
           'v9.1.0 model, old anaerobic script', ...
           'v9.0.2 model, new anaerobic script', ...
           'v9.0.2 model, old anaerobic script'};

fig1 = figure('Name', 'Growth Comparison', 'Position', [100 100 1000 800]);
for k = 1:4
    subplot(2,2,k);
    funcs{k}(models{k});
    title(titles{k});
end
sgtitle('Growth: v9.1.0 vs v9.0.2, new vs old anaerobic script');
saveas(fig1, '..\..\data\testResults\v910_growth.png');

% R2 of growth prediction increased from 0.8256 to 0.9085 as result of all
% curations.

%% Convert to anaerobic
cd('../otherChanges/')
modelAn910 = anaerobicModel(model910);
modelAnOld910 = anaerobicModelOld(model910);
modelAn902 = anaerobicModel(model902);
modelAnOld902 = anaerobicModelOld(model902);

%% Anaerobic flux predictions
cd('../modelTests/');
models = {modelAn910, modelAnOld910, modelAn902, modelAnOld902};
titles  = {'v9.1.0 model, new anaerobic script', ...
           'v9.1.0 model, old anaerobic script', ...
           'v9.0.2 model, new anaerobic script', ...
           'v9.0.2 model, old anaerobic script'};
fig2 = figure('Name', 'Anaerobic Comparison', 'Position', [100 100 1000 800]);
for k = 1:4
    subplot(2,2,k);
    anaerobic_flux_predictions(models{k});
    title(titles{k});
end
sgtitle('Anaerobic flux predictions: v9.1.0 vs v9.0.2, new vs old anaerobic script');
saveas(fig2, '..\..\data\testResults\v910_anaerobic_fluxes.png');

% R2 of flux predictions increased from 0.7858 to 0.9175 as a response to
% all curations

%% Plot
models = {modelAn910, modelAnOld910, modelAn902, modelAnOld902};
titles  = {'v9.1.0 model, new anaerobic script', ...
           'v9.1.0 model, old anaerobic script', ...
           'v9.0.2 model, new anaerobic script', ...
           'v9.0.2 model, old anaerobic script'};
fig3 = figure('Name', 'Anaerobic Comparison', 'Position', [100 100 1000 800]);
exchange = cell(1,4);
for k = 1:4
    subplot(2,2,k);
    exchange{k} = plotAnaerobic(models{k});
    title(titles{k});
end
sgtitle('Anaerobic fluxes (Sjöberg data): v9.1.0 vs v9.0.2, new vs old anaerobic script');
saveas(fig3, '..\..\data\testResults\v910_anaerobic_sjoberg.png');

%% Exchange rate predictions against Sjöberg et al. (2024)
% Ethanol is predicted somewhat above the measured rate; the other three
% products fall within their experimental error.
fprintf('\nExchange rates vs Sjoberg et al. (2024):\n')
for k = 1:4
    fprintf('  %-38s mean rel. error %.4f, %.0f%% within error, NH4+/ATPase %.2f\n', ...
        titles{k}, exchange{k}.meanRelativeError, ...
        100*exchange{k}.fractionWithinError, exchange{k}.ammoniumPerATPase);
end

% The ammonium/ATPase ratio is measured to be close to 1.

%% Pack flux results into table
rxns_reacs=constructEquations(temp_model,temp_model.rxns);
tab=table(temp_model.rxns,temp_model.rxnNames,rxns_reacs,abs(FLUX./v_glc),FLUX,temp_model.grRules);

% [massImbalance, imBalancedMass, imBalancedCharge, imBalancedRxnBool, elements, missingFormulaeBool, balancedMetBool] = checkMassChargeBalance(temp_model);
% 
% tabImbalance = [tab,table(imBalancedRxnBool,imBalancedCharge,imBalancedMass)];
% tabImbalance(~imBalancedRxnBool,:) = [];
% % Filter out exchange and SLIME reactions, these are known to be unbalanced
% tabImbalance(contains(tabImbalance.Var2,'exchange'),:) = [];
% tabImbalance(contains(tabImbalance.Var2,'SLIME'),:) = [];


%% Calculate formula and degree of reduciton of biomass
% [mwRange,metFormulae,elements,metEle]=computeMetFormulae(modelAn,'metMwRange','s_0450','fillMets','none','printLevel',0);
% Biomass_index = find(strcmp(modelAn.metNames,'biomass'));
% %Degree of reduction per element. Order of the elements 'C', 'H', 'O', 'N'
% DR_per_ele = [4, 1, -2, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
% DR_per_Cmol=sum(metEle(Biomass_index,:).*DR_per_ele)/metEle(Biomass_index,1)
% Biomass_formula_Cmol=metEle(Biomass_index,:)/metEle(Biomass_index,1);
% MW_per_Cmol_min = mwRange/metEle(Biomass_index,1)
% 
%

%% Old functions
function R2 = growthOld(model_origin,writeOutput)
% This is for growth test: Fig S4c for yeast8 paper
% here we use several chemostat data: 'N-limited aerobic' 'C-limited
% aerobic' 'C-limited anaerobic' 'N-limited anaerobic'
% when simulating N-limited condition, protein content was rescaled, and
% when simulate anaerobic condtion, heme NADH NADP NADPH NAD were rescaled
% to be 0.

funcDir = fileparts(mfilename('fullpath'));

if nargin<1
    model_origin = loadYeastYaml;
end
if nargin<2
    writeOutput = false;
end

%Load chemostat data:
fid = fopen(fullfile(funcDir,'..','..','data','physiology','chemostatData_Tobias2013.tsv'),'r');
exp_data = textscan(fid,'%f32 %f32 %f32  %f32','Delimiter','\t','HeaderLines',1);
exp_data = [exp_data{1} exp_data{2} exp_data{3} exp_data{4}];
fclose(fid);
exp_data1 = exp_data(1:9,:);
exp_data2 = exp_data(10:20,:);
exp_data3 = exp_data(21:26,:);
exp_data4 = exp_data(27:32,:);

%'N-limited aerobic'
mod_data(1:9,:) = simulateChemostat(model_origin,exp_data(1:9,:),1,'N');
%'C-limited aerobic'
mod_data(10:20,:) = simulateChemostat(model_origin,exp_data(10:20,:),1,'C');
%'C-limited anaerobic'
mod_data(21:26,:) = simulateChemostat(model_origin,exp_data(21:26,:),2,'C');
%'N-limited anaerobic'
mod_data(27:32,:) = simulateChemostat(model_origin,exp_data(27:32,:),2,'N');

% plot the figure
hold on
cols = [215,25,28;253,174,97;171,217,233;44,123,182]/256;
b(1) = plot(exp_data(1:9,4),mod_data(1:9,4),'o','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',cols(2,:));
b(2) = plot(exp_data(10:20,4),mod_data(10:20,4),'s','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',cols(1,:));
b(3) = plot(exp_data(21:26,4),mod_data(21:26,4),'d','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',cols(3,:));
b(4) = plot(exp_data(27:32,4),mod_data(27:32,4),'>','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',cols(4,:));
exp_max = max(exp_data(:,4));
mod_max = max(mod_data(:,4));
lim = max(exp_max,mod_max)+0.05;
xlim([0 lim])
ylim([0 lim])
x=0:0.001:lim;
y = x;
plot(x,y,'--','MarkerSize',6,'Color',[64,64,64]/256)
xlabel('Experimental growth rate [1/h]','FontSize',14,'FontName','Helvetica')
ylabel('In silico growth rate [1/h]','FontSize',14,'FontName','Helvetica')
legend(b,'N-limited aerobic','C-limited aerobic','C-limited anaerobic','N-limited anaerobic','Location','northwest')

% meanerror = sqrt(sum(([exp_data1(:,4);exp_data2(:,4);exp_data3(:,4);exp_data4(:,4)]-[mod_data1(:,4);mod_data2(:,4);mod_data3(:,4);mod_data4(:,4)]).^2)/32)/sqrt(32);
% text(0.25,0.1,['SEM:',num2str(meanerror)])
hold off
R2=corrcoef(exp_data(:,4),mod_data(:,4));
R2=R2(2)^2;
text(0.25,0.1,['R2:',num2str(R2)])

if writeOutput
    saveas(gcf,fullfile(funcDir,'..','..','..','data','testResults','growth.png'));
    fid = fopen(fullfile(funcDir,'..','..','..','data','testResults','growth.md'),'w');
    fprintf(fid,'%s\n','## R2 of growth rate prediction');
    fprintf(fid,'%.4g\n\n',R2);
    fprintf(fid,'%s\n','![Growth curve](growth.png)');
    fclose(fid);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mod_data,solresult] = simulateChemostat(model_origin,exp_data,mode1,mode2)
%Relevant positions:
pos(1) = find(strcmp(model_origin.rxns,'r_1714')); %glc
pos(2) = find(strcmp(model_origin.rxns,'r_1992')); %O2
pos(3) = find(strcmp(model_origin.rxns,'r_1654')); %NH3
pos(4) = find(strcmp(model_origin.rxns,'r_2111')); %growth

%Simulate chemostats:
mod_data = zeros(size(exp_data));
solresult = zeros(length(model_origin.rxns),length(exp_data(:,1)));
if mode1 == 2
    model_origin = anaerobicModelOld(model_origin);
end
if strcmp(mode2,'N')
    % P content in NH3-lim 0.1/h chemostat, per 10.1016/j.femsyr.2005.04.003
    model_origin = scaleBioMassOld(model_origin,'protein',0.28,'',false);
    % Assume that RNA decreased by the same amount (40%)
    model_origin = scaleBioMassOld(model_origin,'RNA',0.0329,'carbohydrate',false);
    model_origin = setParam(model_origin,'ub','r_0472',1000); %Glutamate synthase repressed in excess nitrogen
end
for i = 1:length(exp_data(:,1))
    model_test= model_origin;
    %Fix glucose uptake rate and maximize growth:
    for j = 1:length(exp_data(1,:))-1

        if abs(exp_data(i,j))==1000
            model_test = setParam(model_test,'lb',model_test.rxns(pos(j)),-exp_data(i,j));
        else
            model_test = setParam(model_test,'eq',model_test.rxns(pos(j)),-exp_data(i,j));
        end
    end

    model_test = setParam(model_test,'obj',model_test.rxns(pos(4)),1);
    sol        = solveLP(model_test,1);
    %Store relevant variables:
    try
        mod_data(i,:) = abs(sol.x(pos)');
        solresult(:,i) = sol.x;
    catch
        mod_data(i,:) = 0;
        solresult(:,i) = 0;
    end
end
end

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

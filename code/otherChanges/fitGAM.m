%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = fitGAM(model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = fitGAM(model)

%Load chemostat data:
fid = fopen('../../data/physiology/chemostatData_VanHoek1998.tsv','r');
exp_data = textscan(fid,'%f32 %f32 %f32 %f32','Delimiter','\t','HeaderLines',1);
exp_data = [exp_data{1} exp_data{2} exp_data{3} exp_data{4}];
fclose(fid);

%GAMs to span:
disp('Estimating GAM:')
GAM = 30:5:70;

%1st iteration:
GAM = iteration(model,GAM,exp_data);

%2nd iteration:
GAM = iteration(model,GAM-10:1:GAM+10,exp_data);

%3rd iteration:
GAM = iteration(model,GAM-1:0.1:GAM+1,exp_data);

model = changeGAM(model,GAM);

%Plot fit:
mod_data = simulateChemostat(model,exp_data);
figure
hold on
cols = [0,1,0;0,0,1;1,0,0];
b    = zeros(1,length(exp_data(1,:))-1);
for i = 1:length(exp_data(1,:))-1
    b(i) = plot(mod_data(:,1),mod_data(:,i+1),'Color',cols(i,:),'LineWidth',2);
    plot(exp_data(:,1),exp_data(:,i+1),'o','Color',cols(i,:),'MarkerFaceColor',cols(i,:))
end
title('GAM fitting for growth on glucose minimal media')
xlabel('Dilution rate [1/h]')
ylabel('Exchange fluxes [mmol/gDWh]')
legend(b,'Glucose consumption','O2 consumption','CO2 production','Location','northwest')
hold off

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GAM = iteration(model,GAM,exp_data)

fitting = ones(size(GAM))*1000;

for i = 1:length(GAM)
    %Modify GAM:
    model_i = changeGAM(model,GAM(i));
    
    %Simulate model and calculate fitting:
    mod_data   = simulateChemostat(model_i,exp_data);
    R          = (mod_data - exp_data)./exp_data;
    fitting(i) = sqrt(sum(sum(R.^2)));
    disp(['GAM = ' num2str(GAM(i)) ' -> Error = ' num2str(fitting(i))])
end

%Choose best:
[~,best] = min(fitting);

if best == 1 || best == length(GAM)
    error('GAM found is sub-optimal: please expand GAM search bounds.')
else
    GAM = GAM(best);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = changeGAM(model,GAM)

bioPos = strcmp(model.rxnNames,'biomass pseudoreaction');
for i = 1:length(model.mets)
    S_ix  = model.S(i,bioPos);
    isGAM = sum(strcmp({'ATP [cytoplasm]','ADP [cytoplasm]','H2O [cytoplasm]', ...
        'H+ [cytoplasm]','phosphate [cytoplasm]'},model.metNames{i})) == 1;
    if S_ix ~= 0 && isGAM
        model.S(i,bioPos) = sign(S_ix)*GAM;
    end
end

end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mod_data = simulateChemostat(model,exp_data)

%Relevant positions:
pos(1) = find(strcmp(model.rxnNames,'growth'));
pos(2) = find(strcmp(model.rxnNames,'D-glucose exchange'));
pos(3) = find(strcmp(model.rxnNames,'oxygen exchange'));
pos(4) = find(strcmp(model.rxnNames,'carbon dioxide exchange'));

%Simulate chemostats:
mod_data = zeros(size(exp_data));
for i = 1:length(exp_data(:,1))
    %Fix biomass and minimize glucose:
    model = changeRxnBounds(model,model.rxns(pos(1)),exp_data(i,1),'l');
    model = changeRxnBounds(model,model.rxns(pos(2)),-10,'l');
    model = changeObjective(model,model.rxns(pos(2)),+1);
    sol   = optimizeCbModel(model);
    %Store relevant variables:
    mod_data(i,:) = abs(sol.x(pos)');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

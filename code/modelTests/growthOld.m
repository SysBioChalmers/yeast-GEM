function R2 = growthOld(model_origin,writeOutput)
% This is for growth test: Fig S4c for yeast8 paper
% here we use several chemostat data: 'N-limited aerobic' 'C-limited
% aerobic' 'C-limited anaerobic' 'N-limited anaerobic'
% when simulating N-limited condition, protein content was rescaled, and
% when simulate anaerobic condtion, heme NADH NADP NADPH NAD were rescaled
% to be 0.

funcDir = dbstack('-completenames');
funcDir = regexprep(funcDir(1).file,[funcDir(1).name '\.m'],'');

if nargin<1
    model_origin = loadYeastModel;
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
    saveas(gcf,fullfile(funcDir,'..','..','data','testResults','growth.png'));
    fid = fopen(fullfile(funcDir,'..','..','data','testResults','growth.md'),'w');
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

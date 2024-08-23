function R2=anaerobic_flux_predictions(model)

[vals_flux text_flux]=xlsread('../../data/physiology/anaerobic_flux.xlsx','flux_data');

% %% This is temporary; remove when rxn id is attributed
% text_flux(strcmp(text_flux(:,7),'rxxS'),7)={'r_1239'};

sim_vals=[];


colors = orderedcolors("glow12");

figure;
data_sets=unique(text_flux(:,6));
merged_data=[];
merged_sim=[];
merged_names=[];

for i=1:length(data_sets)

    index_data_set=find(strcmp(text_flux(:,6),data_sets(i)));
    index_glx=strcmp(text_flux(index_data_set,7),'r_1714');
    model = setParam(model,'eq','r_1714',-vals_flux(index_data_set(index_glx),2) );

    %% Solve the LP problem
    res=solveLP(model,1);

    index_model=findRxnIDs(model,text_flux(index_data_set,7));
    exclude_data=(index_model==0);
    index_model(exclude_data)=[];
    index_data_set(exclude_data)=[];

    scaled_sim=abs(-100.*res.x(index_model)./res.x(findRxnIDs(model,'r_1714')));

    data_vals=abs(vals_flux(index_data_set,3));
    merged_data=[merged_data data_vals'];
    merged_sim=[merged_sim scaled_sim'];
    merged_names=[merged_names text_flux(index_data_set,7)'];
    plot(data_vals,scaled_sim,'^','MarkerFaceColor',colors(i,:),'MarkerEdgeColor',colors(i,:));
    hold on;

    % res.x()
end



threshold=30;
R2=corrcoef(merged_sim(merged_data<threshold),merged_data(merged_data<threshold));
R2=R2(2)^2;
mean_relative_error=mean(abs((merged_sim(merged_data<threshold)-merged_data(merged_data<threshold))./merged_data(merged_data<threshold)),'omitnan');

x=0:1:threshold;
y = x;
plot(x,y,'--','MarkerSize',6,'Color',[64,64,64]/256)
ylim([0 threshold])
xlim([0 threshold])
text(5,threshold/2,['mean relative error: ' num2str(mean_relative_error)]);
text(5,threshold/2-10,['R^2: ' num2str(R2)]);

legend(data_sets);
xlabel('Experimental 100 \cdot v_i/v_{Glx}','FontSize',14,'FontName','Helvetica')
ylabel('In silico 100 \cdot v_i/v_{Glx}','FontSize',14,'FontName','Helvetica');

end
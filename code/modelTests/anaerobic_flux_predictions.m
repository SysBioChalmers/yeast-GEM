function flux=anaerobic_flux_predictions(model)
%   flux   struct with the fold-error metrics, the count of measurements
%          that have no comparable ratio, and R2 (kept for comparison
%          with earlier releases; see data/testResults/README.md for why
%          it is not the headline number)

fluxTable = readtable('../../data/physiology/flux_data_anaerobic.tsv','FileType','text');
fluxTable = table2cell(fluxTable);
vals_flux=fluxTable;
text_flux=fluxTable;

sim_vals=[];

colors = orderedcolors("glow12");

data_sets=unique(text_flux(:,6));
merged_data=[];
merged_sim=[];
merged_names=[];

for i=1:length(data_sets)
    % Gather single data set
    idx_data_set = find(strcmp(text_flux(:,6),data_sets(i,1)));
    idx_glc = idx_data_set(strcmp(text_flux(idx_data_set,7),'r_1714'));
    % Set glucose uptake
    model = setParam(model,'eq','r_1714',-cell2mat(vals_flux(idx_glc,4)));
    % Solve the LP problem
    res=solveLP(model,1);
    % Organize data and output
    rxns                        = text_flux(idx_data_set,7);
    include_data                = ismember(rxns,model.rxns);

    rxns(~include_data)         = [];
    idx_data_set(~include_data) = [];
    
    idx_model = getIndexes(model,rxns,'rxns');
    
    scaled_sim=abs(-100.*res.x(idx_model)./res.x(getIndexes(model,'r_1714','rxns')));

    data_vals=abs(cell2mat(vals_flux(idx_data_set,5)));
    merged_data=[merged_data data_vals'];
    merged_sim=[merged_sim scaled_sim'];
    merged_names=[merged_names rxns'];
    plot(data_vals,scaled_sim,'^','MarkerFaceColor',colors(i,:),'MarkerEdgeColor',colors(i,:));
    hold on;

end



threshold=30;
% R2 = coefficient of determination about the line of identity, with the
% experimental values (merged_data) as reference; all data points (no threshold).
R2 = 1 - sum((merged_data-merged_sim).^2)/sum((merged_data-mean(merged_data)).^2);
mean_relative_error = mean(abs((merged_sim-merged_data)./merged_data),'omitnan');

%Fold error. These fluxes span three orders of magnitude, so what matters
%is how far out a prediction is proportionally: 2x out at a flux of 1 is
%the same error as 2x out at a flux of 100, which a difference-based score
%does not say. A ratio needs both values non-zero and pointing the same
%way; pairs that fail either test are counted rather than dropped, since
%predicting no flux where flux was measured, or the opposite direction, is
%a worse error than any finite ratio.
fluxTol = 1e-9;
comparable = abs(merged_data) > fluxTol & abs(merged_sim) > fluxTol & ...
             sign(merged_sim) == sign(merged_data);
logRatio = abs(log10(abs(merged_sim(comparable))./abs(merged_data(comparable))));

flux.medianFoldError = 10^median(logRatio);
flux.meanFoldError   = 10^mean(logRatio);
flux.withinTwoFold   = mean(logRatio < log10(2));
flux.unpredicted     = numel(merged_data) - sum(comparable);
flux.R2              = R2;
flux.meanRelativeError = mean_relative_error;

x=0:1:threshold;
y = x;
plot(x,y,'--','MarkerSize',6,'Color',[64,64,64]/256)
ylim([0 threshold])
xlim([0 threshold])
text(12,threshold/2-10,['mean relative error: ' num2str(mean_relative_error)]);
text(12,threshold/2-5,['R^2: ' num2str(R2)]);

legend(data_sets);
xlabel('Experimental 100 \cdot v_i/v_{Glx}','FontSize',14,'FontName','Helvetica')
ylabel('In silico 100 \cdot v_i/v_{Glx}','FontSize',14,'FontName','Helvetica');

end
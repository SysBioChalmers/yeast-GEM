function plotAnaerobic(modelAn)
modelAn = setParam(modelAn,'eq','r_1714',-23);
res=solveLP(modelAn,1);
FLUX=res.x;

%% Retrieve data for the main products
v_glc=FLUX(getIndexes(modelAn,'r_1714','rxns'),:);
v_eth=FLUX(getIndexes(modelAn,'r_1761','rxns'),:);
v_CO2=FLUX(getIndexes(modelAn,'r_1672','rxns'),:);
v_gly=FLUX(getIndexes(modelAn,'r_1808','rxns'),:);
v_growth=FLUX(getIndexes(modelAn,'r_4041','rxns'),:);
v_AStr = FLUX(getIndexes(modelAn,'r_1115','rxns'));
v_ATPase = FLUX(getIndexes(modelAn,'r_0227','rxns'));
%% Show relative accuracy of main extracellular products
%figure;
%glycerol ethanol Co2
%4.5 ± 0.4  31 ± 2  38 ± 10
data=[4.5 31 38 0.36];
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
hold off
end
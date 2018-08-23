%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = updateMetaboliteFormula(model)
%
% Reads the data file and updates the metabolite formula information in the model
%
% INPUTS:
% model         SBML model used in the simulation
% 
% OUTPUTS:
% model         SBML model used in the simulation (changed)
%
%
% William T. Scott, Jr.
% Last Update: 2018-08-23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function model = updateMetaboliteFormula(~)
% Load model
cd ..
model = loadYeastModel;
 
%Correct some of the metabolite fomulae from most recent info
%Load data:
fid = fopen('../../ComplementaryData/modelCuration/Missingmetaboliteformulas.tsv','r');
metaboliteData = textscan(fid,'%s %s %s %s %f32 %s','Delimiter','\t','HeaderLines',1);
fclose(fid);
 
for i = 1:length(metaboliteData{1})
    for j = 1:length(model.mets)
        %Find name:
        if startsWith(model.metNames{j},[metaboliteData{2}{i} ' ['])    %metabolite name
            metName = model.metNames{j};
            comp    = metName(strfind(metName,' ['):end);
            metName = [metaboliteData{2}{i} comp];
            model.metNames{j} = metName;
        end
        
        %Update other fields:
        if startsWith(model.metNames{j},[metaboliteData{2}{i} ' ['])    %new name
            model.model.metFormulas{j} = metaboliteData{3}{i};  %new formulas    
        end
    end
end
 
% Save model
cd ..
saveYeastModel(model)
cd modelCuration
 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
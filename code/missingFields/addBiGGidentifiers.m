function addBiGGidentifiers
% addBiGGidentifiers
%   Loads the model, adds all BiGG identifiers for both metabolites and
%   reactions, and saves the model back.
%

% Load model:
cd ..
initCobraToolbox
model = loadYeastModel;

% Add metabolite ids:
[~,met_dic] = xlsread('../data/databases/BiGGmetDictionary.csv');
model.metBiGGID = cell(size(model.mets));
for i = 1:size(met_dic,1)
    model.metBiGGID(strcmp(model.mets,met_dic{i,1})) = met_dic(i,2);
end

% Add reaction ids:
[~,rxn_dic] = xlsread('../data/databases/BiGGrxnDictionary.csv');
model.rxnBiGGID = cell(size(model.rxns));
for i = 1:size(rxn_dic,1)
    model.rxnBiGGID(strcmp(model.rxns,rxn_dic{i,1})) = rxn_dic(i,2);
end

% Save model:
saveYeastModel(model);
cd missingFields

end

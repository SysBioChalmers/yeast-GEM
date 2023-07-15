
% Load Model:
cd ..
model = loadYeastModel;

% Migrate KEGG pathways from 'subSystems' to 'rxnKEGGpathways'
model.rxnKEGGPathways = cell(size(model.rxns));
for i = 1:length(model.rxns)
    KEGGPathways = model.subSystems{i};
    for j = 1:length(KEGGPathways)
        if ~isempty(KEGGPathways{j})
            KEGGPathways{j} = KEGGPathways{j}(1:8);
        end
    end
    model.rxnKEGGPathways{i} = strjoin(KEGGPathways,';');
end
model.subSystems = cell(size(model.rxns));

% Save Model:
saveYeastModel(model);
cd otherChanges

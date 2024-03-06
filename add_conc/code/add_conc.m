% load model and YMDB data
modelWithConc = importModel('../../model/yeast-GEM.xml');
data = readtable('../data/allConcData.xlsx');
data = table2cell(data);
metConc.allConc = cell(length(modelWithConc.metNames), 1);
metConc.maxConc = cell(length(modelWithConc.metNames), 1);
metConc.minConc = cell(length(modelWithConc.metNames), 1);

met.name = modelWithConc.metNames;
met.kegg = {};
for i = 1:length(met.name)
    idx = find(strcmp('kegg.compound',modelWithConc.metMiriams{i, 1}.name));
    if idx
        met.kegg{i, 1} = modelWithConc.metMiriams{i, 1}.value{idx, 1};
    end
    clear idx
end

% add conc into model
metWithConc = {};
for j = 2:870
    keggId = data{2, j};
    idx = find(strcmp(keggId, met.kegg));
    if idx
        metWithConc{end+1} = keggId;
        for id= 1:length(idx)
            metConc.allConc{idx(id), 1} = data{3, j};
            metConc.maxConc{idx(id), 1} = data{4, j};
            metConc.minConc{idx(id), 1} = data{5, j};
        end
    end
end
fprintf('In total, %d mets have concentration range.', length(metWithConc))
modelWithConc.metConc = metConc;
save('../modelWithConc.mat', 'modelWithConc');
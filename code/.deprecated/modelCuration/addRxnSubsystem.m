% This Function is for adding new subsystem into model.
% Input: model, keggPathways.tsv.
% NOTE: Before run the codes below, the file should be manually editted.
%       COBRA required.
%       New subsystem should be in .tsv format.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load model
cd ..
model = loadYeastModel;

%newsubsystem:
fid = fopen('../data/modelCuration/v8_6_0/keggPathways.tsv');
sub_data = textscan(fid,'%s %s %s','Delimiter','\t','HeaderLines',1);
id = sub_data{1};
subSystems = sub_data{2};
fclose(fid);
for i = 1:length(model.rxns)
    if i <= length(id)
        subIndex = find(strcmp(model.rxns, id(i)));
        model.subSystems{subIndex, 1} = subSystems{i};
    else
        model.subSystems{i, 1} = 'nan';
    end
end

% Save model:
cd ..
saveYeastModel(model)
cd modelCuration

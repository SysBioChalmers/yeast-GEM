cd ../code/
model = loadYeastModel;
cd ../add_uniprot_id/
fid = fopen('SGD_with_Uniprot.csv');
uniprot = textscan(fid,'%q %q %q','Delimiter',',','HeaderLines',1);
for i = 1:length(model.genes)
    ids.name{1, 1} = char("Uniprot ID");
    ind = find(strcmp(model.genes{i, 1}, uniprot{1, 1}));
    ids.value{1, 1} = uniprot{1, 3}{ind, 1};
    model.geneMiriams{i, 1} = ids;
end
cd ../code/
saveYeastModel(model)
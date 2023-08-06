cd ../code/
model = loadYeastModel;
cd ../add_uniprot_id/
fid = fopen('SGD_with_Uniprot.csv');
uniprot = textscan(fid,'%q %q %q','Delimiter',',','HeaderLines',1);
for i = 1:length(model.genes)
    gM = model.geneMiriams(i);
    gM = gM{1,1};
    if ~isa(gM, 'struct')
        ids.name{1, 1} = char('kegg.genes');
        ids.name{2, 1} = char('ncbigene');
        ids.name{3, 1} = char('refseq');
        ids.name{4, 1} = char('uniprot');
        ids.name{5, 1} = char('ncbiprotein');
        ind = find(strcmp(model.genes{i, 1}, uniprot{1, 1}));
        ids.value{4, 1} = uniprot{1, 3}{ind, 1};
        ids.value{1, 1} = char('NaN');
        ids.value{2, 1} = char('NaN');
        ids.value{3, 1} = char('NaN');
        ids.value{5, 1} = char('NaN');
        model.geneMiriams{i, 1} = ids;
        clear ids
    end
end
cd ../code/
saveYeastModel(model)
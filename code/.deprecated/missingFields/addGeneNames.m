%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = addGeneNames(model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = addGeneNames(model)

%Correct gene ids for proper matching:
model.genes = strrep(model.genes,'_','-');

%Load data:

fid  = fopen('../../data/databases/swissprot.tsv','r');
data = textscan(fid,'%s %s %s %s %f32 %s','Delimiter','\t','HeaderLines',1);
fclose(fid);
swiss.genes = data{3};
for i = 1:length(swiss.genes)
    swiss.genes{i} = strsplit(swiss.genes{i},' ');
end

%Get gene names:
for i = 1:length(model.genes)
    for j = 1:length(swiss.genes)
        if sum(strcmp(swiss.genes{j},model.genes{i})) > 0
            model.geneNames{i} = swiss.genes{j}{1}; %First occurence is the gene name
        end
    end
    disp(['Adding gene names: Ready with gene #' int2str(i)])
end

%Save model:
cd ..
saveYeastModel(model)
cd missingFields

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
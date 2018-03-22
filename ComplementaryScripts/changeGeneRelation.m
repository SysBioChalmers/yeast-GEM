%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correct gene relation (and add some new gene from isce926)
% changeGeneAssociation.m is a function from cobra
% March 22, 2018 by Hongzhong
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = readCbModel('yeastGEM.xml');
% correct some gene relations based on isce926
%fid = fopen('../../ComplementaryData/correct_gene_relation_isce926.tsv');
fid = fopen('correct_gene_relation_isce926.tsv');
correct_gene_relation = textscan(fid,'%s %s %s %s %s','Delimiter','\t','HeaderLines',1);
fclose(fid);

ss1 = length(correct_gene_relation{1})
oldGPRrxns = zeros(ss1,1)
for  i = 1 : ss1
    oldGPRrxns(i) = find(strcmp(model.rxns, correct_gene_relation{1}{i}));%Find all reactions that have the old GPR
    model = changeGeneAssociation(model, model.rxns{oldGPRrxns(i)}, correct_gene_relation{4}{i});
end

% add new genes based on isce926
%fid1 = fopen('../../ComplementaryData/newGene_from_isce926.tsv');
fid1 = fopen('newGene_from_isce926.tsv');
newGene_from_isce926 = textscan(fid1,'%s %s %s %s %s','Delimiter','\t','HeaderLines',1);
fclose(fid1);


ss2 = length(newGene_from_isce926{1})
oldGPRrxns = zeros(ss2,1)
for  i = 1 : ss2
    oldGPRrxns(i) = find(strcmp(model.rxns, newGene_from_isce926{1}{i}));%Find all reactions that have the old GPR
    model = changeGeneAssociation(model, model.rxns{oldGPRrxns(i)}, newGene_from_isce926{4}{i});
end


% add gene standard name for new gene from isce926
%fid2 = fopen('../../ComplementaryData/yeast_gene_annotation_SGD.tsv');
fid2 = fopen('yeast_gene_annotation_SGD.tsv');
yeast_gene_annotation = textscan(fid2,'%s %s','Delimiter','\t','HeaderLines',1);
fclose(fid2);

ss3 = length(model.genes)
genePosition = zeros(ss3,1)
for i = 1: ss3
    if ~isempty(find(strcmp(yeast_gene_annotation{1}, model.genes{i})))
        genePosition(i) = find(strcmp(yeast_gene_annotation{1}, model.genes{i}))
        model.geneNames{i} = yeast_gene_annotation{2}{genePosition(i)}
    else
        genePosition(i) = nan
    end
end

% add protein name for new gene from isce926
ss4 = length(model.genes)
proteinName = cell(ss4,1)
for i = 1:ss4
proteinName{i} = strcat('COBRARProtein',num2str(i))
model.proteins{i} = proteinName{i}
end


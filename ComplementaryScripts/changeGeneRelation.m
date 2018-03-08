%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correct gene relation (and add some new gene from isce926)
% changeGeneAssociation.m is a function from cobra
% Feb 7, 2018 by Hongzhong
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model2 = readCbModel('yeastGEM.xml');
model = model2; %import the latest model

GPRsReplace = correctgenerelationisce926 % import as string
GPRsReplace = cellstr(GPRsReplace) % change into cell format

newGene = newgeneisce926 % import as string
newGene =  cellstr(newGene) % change into cell format

gene_name = yeastgeneannotationSGD % import as string

gene_name(:,2) % gene systematic name in SGD
gene_standard_name = gene_name(:,4) % gene standard name
% correct some gene relation based on isce926
for  i = 1 : size(GPRsReplace, 1)
    oldGPRrxns = find(strcmp(model.rules, GPRsReplace{i, 1}));%Find all reactions that have the old GPR
    for j = 1:length(oldGPRrxns)
        model = changeGeneAssociation(model, model.rxns{oldGPRrxns(j)}, GPRsReplace{i, 2});
    end
end



% add some new gene based on isce926
t1 = length(newGene)
ss = zeros(t1,1)
for  i = 1:t1
    ss(i) = strmatch(newGene(i,1), model.rxns)
    model = changeGeneAssociation(model,model.rxns{ss(i)}, newGene{i,2})
end


% add gene standard name for gene
t2 = length(model.genes)
ss1 = zeros(t2,1)
for i = 1: t2
    if ~isempty(find(strcmp(gene_name(:,2), model.genes{i})))
        ss1(i) = find(strcmp(gene_name(:,2), model.genes{i}))
        model.geneNames{i} = gene_standard_name{ss1(i)}
    else
        ss1(i) = nan
    end
end

% add protein name for gene
t2 = length(model.genes)
ss2 = cell(t2,1)
for i = 1:t2
ss2{i} = strcat('COBRAProtein',num2str(i))
model.proteins{i} = ss2{i}
end

saveYeastModel(model)

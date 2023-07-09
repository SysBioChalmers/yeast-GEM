%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addiSce926changes
% Add changes from the model iSce926 + manual curation on those changes
% iSce926 source: http://www.maranasgroup.com/submission_models/iSce926.htm
%
% NOTE: changeGeneAssociation.m is a function from cobra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load model:
model = readCbModel('../../model/yeast-GEM.xml');

%Correct some gene relations based on isce926:
fid      = fopen('../../data/modelCuration/iSce926curatedGeneRules.tsv');
newRules = textscan(fid,'%s %s %s %s %s','Delimiter','\t','HeaderLines',1);
fclose(fid);

ss1 = length(newRules{1});
oldGPRrxns = zeros(ss1,1);
for i = 1:ss1
    %Find all reactions that have the old GPR:
    oldGPRrxns(i) = find(strcmp(model.rxns, newRules{1}{i}));
    model         = changeGeneAssociation(model, model.rxns{oldGPRrxns(i)}, newRules{4}{i});
end

%Add new genes based on isce926:
fid1     = fopen('../../data/modelCuration/iSce926newGenes.tsv');
newGenes = textscan(fid1,'%s %s %s %s %s','Delimiter','\t','HeaderLines',1);
fclose(fid1);

ss2 = length(newGenes{1});
oldGPRrxns = zeros(ss2,1);
for i = 1:ss2
    %Find all reactions that have the old GPR:
    oldGPRrxns(i) = find(strcmp(model.rxns, newGenes{1}{i}));
    model         = changeGeneAssociation(model, model.rxns{oldGPRrxns(i)}, newGenes{4}{i});
end

%Add gene standard name for new gene from isce926:
fid2 = fopen('../../data/databases/SGDgeneNames.tsv');
SGD  = textscan(fid2,'%s %s','Delimiter','\t','HeaderLines',1);
fclose(fid2);

ss3 = length(model.genes);
genePosition = NaN(ss3,1);
model.geneNames = model.genes;
for i = 1:ss3
    if ismember(model.genes{i},SGD{1})
        genePosition(i) = find(strcmp(SGD{1}, model.genes{i}));
        if ~isempty(SGD{2}{genePosition(i)})
            model.geneNames{i} = SGD{2}{genePosition(i)};
        end
    end
end

%Add protein name for new gene from isce926:
ss4         = length(model.genes);
proteinName = cell(ss4,1);

for i = 1:ss4
    proteinName{i}    = strcat('COBRAProtein',num2str(i));
    model.proteins{i} = proteinName{i};
end

cd ..
saveYeastModel(model)
cd modelCuration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
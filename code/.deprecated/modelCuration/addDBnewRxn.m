% This Function is for adding new database annotated metabolites/reactions into model.
% Add changes from the database new anootation for new reactions and new metabolites and new genes related + manual curation on those changes
% Input: model, DBnewrRxnMatrix.tsv,DBnewRxnMetAnnotation.tsv,SGDgeneNames.tsv.
% As for the reference of new genes and reactions related, please find detailed information in the /data/databases/DBnewGeneAnnotation.tsv
% NOTE: changeGeneAssociation.m is a function from cobra
%       Extract model info from .tsv format.
%       Before run the codes below, the file should be manually editted.
%       COBRA required.
%       New reaction should be in .tsv format.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load model
cd ..
model = loadYeastModel;

%newreaction:
fid = fopen('../data/modelCuration/DBnewrRxnMatrix.tsv');
newreaction    = textscan(fid,'%s %s %s %s %s','Delimiter','\t','HeaderLines',1);
matrix.rxnIDs  = newreaction{1};
matrix.metcoef = cellfun(@str2num, newreaction{2});
matrix.metIDs  = newreaction{3};
matrix.mettype = newreaction{4};
matrix.metcompartments = newreaction{5};
fclose(fid);

%change coefficient
for i=1:length(matrix.rxnIDs)
    if strcmp(matrix.mettype(i),'reactant')
        matrix.metcoef(i) = matrix.metcoef(i)*-1;
        %matrix.metcoef_temp(i) = matrix.metcoef_temp(i)*-1
    end
end

%change compartments
CONValldata = cat(2,model.compNames,model.comps);
lbracket    = ' [' ;%  space
llbracket   = '[';
rbrackets   = ']';
space       = ' ';
[m, n]      = size(CONValldata);
for i = 1:m
    aa = CONValldata(i,1);
    aa = char(aa);
    for j=1:length(matrix.rxnIDs)
        bb = matrix.metcompartments(j,1);
        bb = char(bb);
        if strcmp(bb,aa)
            matrix.Newcomps(j,1) = CONValldata(i,2);
        end
    end
end
for i=1:length(matrix.rxnIDs)
    matrix.metnames(i) = strcat(matrix.metIDs(i),lbracket,matrix.metcompartments(i),rbrackets);
    matrix.Newcomps(i) = strcat(llbracket,matrix.Newcomps(i),rbrackets);
end

%mapping mets to model.metnames, get s_ index for new mets
cd modelCuration/
for j = 1:length(matrix.metnames)
    [~,metindex] = ismember(matrix.metnames(j),model.metNames);
    if metindex ~= 0
        matrix.mets(j) = model.mets(metindex);
    elseif metindex == 0
        newID = getNewIndex(model.mets);
        matrix.mets(j) = strcat('s_',newID,matrix.Newcomps(j));
        model = addMetabolite(model,char(matrix.mets(j)), ...
                              'metName',matrix.metnames(j));
    end
end

% add met annotation
fid = fopen('../../data/modelCuration/DBnewRxnMetAnnotation.tsv');
newmet_annot = textscan(fid,'%s %s %s %s %s %s %s','Delimiter','\t','HeaderLines',1);
newmet.metNames    = newmet_annot{1};
newmet.metFormulas = newmet_annot{2};
newmet.metCharges  = cellfun(@str2num, newmet_annot{3});
newmet.metKEGGID   = newmet_annot{5};
newmet.metChEBIID  = newmet_annot{6};
newmet.metNotes    = newmet_annot{7};
fclose(fid);

for i = 1:length(newmet.metNames)
    [~,metID] = ismember(newmet.metNames(i),model.metNames);
    if metID ~= 0
        model.metFormulas{metID} = newmet.metFormulas{i};
        model.metCharges(metID)  =  newmet.metCharges(i);
        model.metKEGGID{metID}   = newmet.metKEGGID{i};
        model.metChEBIID{metID}  = newmet.metChEBIID{i};
        model.metNotes{metID}    = ['NOTES: added after new annotation (PR #142); ',newmet.metNotes{i}];
    end
end

%Load rxnProp(rev and GPR)
fid = fopen('../../data/modelCuration/DBnewRxnProp.tsv');
rev = textscan(fid,'%s %s %s %s %s %s %s','Delimiter','\t','HeaderLines',1);
newrxn.ID  = rev{1};
newrxn.Rev = cellfun(@str2num, rev{2});
newrxn.GPR = rev{3};
newrxn.rxnNames     = rev{4};
newrxn.rxnECNumbers = rev{5};
newrxn.subSystems   = rev{6};
newrxn.rxnKEGGID    = rev{7};
fclose(fid);

%add new reactions according to rev ID. Met Coef need to be in the column,
%not a row. Coef should be double, which was converted at the import
%section.
EnergyResults     = {};
MassChargeresults = {};
RedoxResults      = {};
for i = 1:length(newrxn.ID)
    cd ../otherchanges
    newID   = getNewIndex(model.rxns);
    cd ../modelCuration
    j = find(strcmp(matrix.rxnIDs,newrxn.ID{i}));
    Met = matrix.mets(j);
    Coef = transpose(matrix.metcoef(j));
    [model,rxnIndex] = addReaction(model,...
                        ['r_' newID],...
                        'reactionName', newrxn.ID{i},...
                        'metaboliteList',Met,...
                        'stoichCoeffList',Coef,...
                        'reversible',newrxn.Rev(i,1),...
                        'geneRule',newrxn.GPR{i},...
                        'checkDuplicate',1);
    [EnergyResults,RedoxResults] = CheckEnergyProduction(model,{['r_' newID]},EnergyResults,RedoxResults);
    [MassChargeresults] = CheckBalanceforSce(model,{['r_' newID]},MassChargeresults);
    if isempty(rxnIndex)
        rxnIndex = strcmp(model.rxns,['r_' newID]);
    end
    % Add rxn annotation:
    model.rxnNames{rxnIndex}      = newrxn.rxnNames{i};
    model.rxnECNumbers(rxnIndex)  = newrxn.rxnECNumbers(i);
    model.rxnKEGGID(rxnIndex)     = newrxn.rxnKEGGID(i);
    model.rxnConfidenceScores(rxnIndex) = 2;   %reactions added 
    model.rxnNotes{rxnIndex} = 'NOTES: added after new annotation (PR #142)';
end


% add gene standard name for new genes
fid = fopen('../../data/databases/SGDgeneNames.tsv');
yeast_gene_annotation = textscan(fid,'%s %s','Delimiter','\t','HeaderLines',1);
fclose(fid);

geneIndex = zeros(1,1);
for i = 1: length(model.genes)
    geneIndex = strcmp(yeast_gene_annotation{1}, model.genes{i});
    if sum(geneIndex) == 1 && ~isempty(yeast_gene_annotation{2}{geneIndex})
        model.geneNames{i} = yeast_gene_annotation{2}{geneIndex};
    else
        model.geneNames{i} = model.genes{i};
    end
end

% Add protein name for genes
for i = 1:length(model.genes)
    model.proteins{i} = strcat('COBRAProtein',num2str(i));
end

% Save model:
cd ..
saveYeastModel(model)
cd modelCuration

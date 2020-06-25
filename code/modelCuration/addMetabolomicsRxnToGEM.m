% This Function is for adding Metabolomics related metabolites/reactions into model.
% Input: model, Metabolomics_newRxnMatrix.tsv,Metabolomics_newRxnProp.tsv,Metabolomics_newRxnMetAnnotation.tsv.
%       Extract model info from .tsv format.
%       Before run the codes below, the file should be manually editted.
%       COBRA required.
%       New reaction should be in .tsv format.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%newreaction:
% Load model
cd ..
model = loadYeastModel;
%load rxn matrix
fid = fopen('../data/modelCuration/Metabolomics_newRxnMatrix.tsv');
newreaction = textscan(fid,'%s %s %s %s %s','Delimiter','\t','HeaderLines',1);
matrix.rxnIDs    = newreaction{1};
matrix.metcoef    = cellfun(@str2num, newreaction{2});
matrix.metIDs = newreaction{3};
matrix.mettype = newreaction{4};
matrix.metcompartments = newreaction{5};
fclose(fid);
%load rxn prop, which include rxn rev, rxn names, rxn ec numbers, rxn
%KEGGID and rxn notes.
fid  = fopen('../data/modelCuration/Metabolomics_newRxnProp.tsv','r');
rev = textscan(fid,'%s %s %s %s %s %s %s','Delimiter','\t','HeaderLines',1);
newrxn.ID  = rev{1};
newrxn.Rev = cellfun(@str2num, rev{2});
newrxn.rxnNames = rev{3};
newrxn.rxnECNumbers = rev{4};
newrxn.rxnKEGGID = rev{5};
newrxn.rxnNotes  = rev{6};
newrxn.rxnMetaNetXID = newrxn.ID;
for i = 1:length(newrxn.rxnMetaNetXID)
    if ~startsWith(newrxn.rxnMetaNetXID{i},'MNXR')
        newrxn.rxnMetaNetXID{i} = '';
    end
end

fclose(fid);

% Change coefficients for reactants:
for i=1:length(matrix.rxnIDs)
    if strcmp(matrix.mettype(i),'reactant')
        matrix.metcoef(i) = matrix.metcoef(i)*-1;
    end
end

% Change compartments:
CONValldata = cat(2,model.compNames,model.comps);
lbracket    = ' [' ;
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

% Map mets to model.metnames, get s_index for new mets:
cd otherChanges
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


% Add metabolite data:
fid = fopen('../../data/modelCuration/Metabolomics_newRxnMetAnnotation.tsv');
newmet_annot         = textscan(fid,'%s %s %s %s %s %s %s %s','Delimiter','\t','HeaderLines',1);
newmet.metNames      = newmet_annot{1};
newmet.metFormulas   = newmet_annot{2};
newmet.metCharges    = cellfun(@str2num, newmet_annot{3});
newmet.metKEGGID     = newmet_annot{5};
newmet.metChEBIID    = newmet_annot{6};
newmet.metMetaNetXID = newmet_annot{7};


fclose(fid);
for i = 1:length(newmet.metNames)
    [~,metID] = ismember(newmet.metNames(i),model.metNames);
    if metID ~= 0
        model.metFormulas{metID}   = newmet.metFormulas{i};
        model.metCharges(metID)    = newmet.metCharges(i);
        model.metKEGGID{metID}     = newmet.metKEGGID{i};
        model.metChEBIID{metID}    = newmet.metChEBIID{i};
        model.metMetaNetXID{metID} = newmet.metMetaNetXID{i};
        model.metNotes{metID}      = 'added from metabolomics data (PR #156)';
    end
end

% Add new reactions according to rev ID: Met Coef needs to be a column, not
% a row. Coef should be a double, which was converted at the import section
EnergyResults     = {};
MassChargeresults = {};
RedoxResults      = {};
if ~isfield(model,'rxnMetaNetXID')
    model.rxnMetaNetXID = cell(size(model.rxns));
end
for i = 1:length(newrxn.ID)
    cd ../otherChanges
    newID = getNewIndex(model.rxns);
    j     = find(strcmp(matrix.rxnIDs,newrxn.ID{i}));
    Met   = matrix.mets(j);
    Coef  = transpose(matrix.metcoef(j));
    [model,rxnIndex] = addReaction(model, ['r_' newID],...
        'reactionName', newrxn.ID{i},...
        'metaboliteList',Met,...
        'stoichCoeffList',Coef,...
        'reversible',newrxn.Rev(i,1),...
        'checkDuplicate',1);
    cd ../modelCuration
    [EnergyResults,RedoxResults] = CheckEnergyProduction(model,{['r_' newID]},EnergyResults,RedoxResults);
    [MassChargeresults] = CheckBalanceforSce(model,{['r_' newID]},MassChargeresults);
    if isempty(rxnIndex)
        rxnIndex = strcmp(model.rxns,['r_' newID]);
    end
    % Add rxn annotation:
    model.rxnNames{rxnIndex}      = newrxn.rxnNames{i};
    model.rxnECNumbers(rxnIndex)  = newrxn.rxnECNumbers(i);
    model.rxnKEGGID(rxnIndex)     = newrxn.rxnKEGGID(i);
    model.rxnMetaNetXID(rxnIndex) = newrxn.rxnMetaNetXID(i);
    model.rxnConfidenceScores(rxnIndex) = 0;   %reactions added for metabolomics data
    model.rxnNotes{rxnIndex} = 'metabolites observed in metabolomics data (PR #156)';
end


% Save model:
cd ..
saveYeastModel(model)
cd modelCuration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


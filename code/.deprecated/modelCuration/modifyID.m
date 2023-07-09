% This script modifies reaction and metabolite IDs based on data from modifyID.tsv
% 
% modifyID.tsv includes details on changes to current rxn/metIDs, addition of new
% rxn/metIDs or notation of alternative rxn/metIDs after manual curation of
% deprecated IDs and selected unbalanced reactions in the model
%
% Inputs: model and modify.tsv
%

%Load model
cd ..
model = loadYeastModel;

%Load modifyID.tsv files
fid = fopen('../data/modelcuration/modifyID.tsv');
format = repmat('%s ',1,16);
format = strtrim(format);
temp = textscan(fid,format,'Delimiter','\t','HeaderLines',0);
for i = 1:length(temp)
    curationfile(:,i) = temp{i}; %use {} instead of () for cell array
end
commentLines = startsWith(curationfile(:,1),'#');
curationfile(commentLines,:) = [];
fclose(fid);

%Different data sets are 'bracketed' i.e. separated by % - use for indexing
brackets = startsWith(curationfile(:,1),'%');
idx = find(brackets);
idx_rxnID = idx(1)+1:idx(2)-1;
idx_metID = idx(2)+1:idx(3)-1;

%Separate data into various cell arrays
modifyrxnID = curationfile(idx_rxnID,1:10);
modifymetID = curationfile(idx_metID,1:16);

%% Update rxnMetaNetXID

rxn = modifyrxnID(:,1);
currentID = modifyrxnID(:,2);
newID = modifyrxnID(:,3);
alternativeID = modifyrxnID(:,4);

for i = 1:length(currentID)
    idx_rxn = find(ismember(model.rxns,rxn(i)));
    %Check for new ID, if present = replace current ID/add ID to model
    %if blank = remove current ID in model
    if ~ismember(newID(i),'[]') && contains(newID(i),'MNXR')
        model.rxnMetaNetXID(idx_rxn) = newID(i);
        model.rxnNotes(idx_rxn) = join([model.rxnNotes(idx_rxn),'| MetaNetX ID curated (PR #220)']);
    elseif ismember(newID(i),'Blank')
        model.rxnMetaNetXID(idx_rxn) = {''};
        model.rxnNotes(idx_rxn) = join([model.rxnNotes(idx_rxn),'| MetaNetX ID curated (PR #220)']);
    elseif ~ismember(newID(i),'[]') && ~contains(newID(i),'MNXR')
        warning('Check for error in %s under rxnID curation data of the tsv file', string(rxn(i)));
    end
    
    %Check for alternative ID, if present = add to rxnNotes
    if ~ismember(alternativeID(i),'[]') && contains(alternativeID(i),'MNXR')
        model.rxnNotes(idx_rxn) = join([model.rxnNotes(idx_rxn),'| alternative MetaNetX ID',alternativeID(i),'(PR #220)']);
    elseif ~ismember(alternativeID(i),'[]') && ~contains(alternativeID(i),'MNXR')
        warning('Check for error in %s under rxnID curation data of the tsv file', string(rxn(i)));
    end
end
        
%% Update rxnKEGGID

currentID = modifyrxnID(:,5);
newID = modifyrxnID(:,6);
alternativeID = modifyrxnID(:,7);

for i = 1:length(currentID)
    idx_rxn = find(ismember(model.rxns,rxn(i)));
    %Check for new ID, if present = replace current ID/add ID to model
    %if blank = remove current ID in model
    if ~ismember(newID(i),'[]') && contains(newID(i),'R')
        model.rxnKEGGID(idx_rxn) = newID(i);
        model.rxnNotes(idx_rxn) = join([model.rxnNotes(idx_rxn),'| KEGG ID curated (PR #220)']);
    elseif ismember(newID(i),'Blank')
        model.rxnKEGGID(idx_rxn) = {''};
        model.rxnNotes(idx_rxn) = join([model.rxnNotes(idx_rxn),'| KEGG ID curated (PR #220)']);
    elseif ~ismember(newID(i),'[]') && ~contains(newID(i),'R')
        warning('Check for error in %s under rxnID curation data of the tsv file', string(rxn(i)));
    end
    
    %Check for alternative ID, if present = add to rxnNotes
    if ~ismember(alternativeID(i),'[]') && contains(alternativeID(i),'R')
        model.rxnNotes(idx_rxn) = join([model.rxnNotes(idx_rxn),'| alternative KEGG ID',alternativeID(i),'(PR #220)']);
    elseif ~ismember(alternativeID(i),'[]') && ~contains(alternativeID(i),'R')
        warning('Check for error in %s under rxnID curation data of the tsv file', string(rxn(i)));
    end
end

%remove whitespace(s) when adding notes into rxnNotes
model.rxnNotes(:) = strtrim(model.rxnNotes(:));

%additional changes to rxnNames in r_2117, r_4254 and r_4255
[~,idx] = ismember('r_4254',model.rxns); %rxnFormula: NADH [cytoplasm] + 2 oxygen [cytoplasm] + 2 nitric oxide [cytoplasm]  -> H+ [cytoplasm] + NAD [cytoplasm] + 2 nitrate [cytoplasm]
model.rxnNames(idx) = {'nitric oxide, NADH2:oxygen oxidoreductase'};
[~,idx] = ismember('r_4255',model.rxns); %rxnFormula: NADPH [cytoplasm] + 2 oxygen [cytoplasm] + 2 nitric oxide [cytoplasm]  -> H+ [cytoplasm] + NADP(+) [cytoplasm] + 2 nitrate [cytoplasm]
model.rxnNames(idx) = {'nitric oxide, NADPH2:oxygen oxidoreductase'};
[~,idx] = ismember('r_2117',model.rxns); %rxnFormula: L-phenylalanine [cytoplasm] + pyruvate [cytoplasm]  <=> keto-phenylpyruvate [cytoplasm] + L-alanine [cytoplasm]
model.rxnECNumbers(idx) = {'2.6.1.58; 2.6.1.7'};

%% Update metMetaNetXID

met = modifymetID(:,1);
[~,idx_met] = ismember(met,model.mets);
modelR = ravenCobraWrapper(model);
metNames = modelR.metNames(idx_met);
currentID = modifymetID(:,2);
newID = modifymetID(:,3);
alternativeID = modifymetID(:,4);

for i = 1:length(metNames)
    idx_met = find(ismember(modelR.metNames,metNames(i)));
    %Check for new ID, replace current ID/add new ID into model
    %if blank = remove current ID in model
    for j = 1:length(idx_met)
        if ~ismember(newID(i),'[]') && contains(newID(i),'MNXM')
            model.metMetaNetXID(idx_met(j)) = newID(i);
        model.metNotes(idx_met(j)) = join([model.metNotes(idx_met(j)),'| MetaNetX ID curated (PR #220)']);
        elseif ismember(newID(i),'Blank')
            model.metMetaNetXID(idx_met(j)) = {''};
        model.metNotes(idx_met(j)) = join([model.metNotes(idx_met(j)),'| MetaNetX ID curated (PR #220)']);
        elseif ~ismember(newID(i),'[]') && ~contains(newID(i),'MNXM')
            warning('Check for error in %s under metID curation data of the tsv file', string(met(i)));
        end
        
        if ~ismember(alternativeID(i),'[]') && contains(alternativeID(i),'MNXM')
            model.metNotes(idx_met(j)) = join([model.metNotes(idx_met(j)),'| alternative MetaNetX ID',alternativeID(i),'(PR #220)']);
        elseif ~ismember(alternativeID(i),'[]') && ~contains(alternativeID(i),'MNXM')
            warning('Check for error in %s under metID curation data of the tsv file', string(met(i)));
        end
    end
end 

%% Update metKEGGID

currentID = modifymetID(:,5);
newID = modifymetID(:,6);
alternativeID = modifymetID(:,7);

for i = 1:length(metNames)
    idx_met = find(ismember(modelR.metNames,metNames(i)));
    %Check for new ID, replace current ID/add new ID into model
    %if blank = remove current ID in model
    for j = 1:length(idx_met)
        if ~ismember(newID(i),'[]') && contains(newID(i),'C')
            model.metKEGGID(idx_met(j)) = newID(i);
            model.metNotes(idx_met(j)) = join([model.metNotes(idx_met(j)),'| KEGG ID curated (PR #220)']);
        elseif ismember(newID(i),'Blank')
            model.metKEGGID(idx_met(j)) = {''};
            model.metNotes(idx_met(j)) = join([model.metNotes(idx_met(j)),'| KEGG ID curated (PR #220)']);
        elseif ~ismember(newID(i),'[]') && contains(newID(i),'G')
            warning('new KEGGID %s not added as it contains G and does not fulfil SBML format', string(newID(i)));
        elseif ~ismember(newID(i),'[]') && ~contains(newID(i),'C')
            warning('Check for error in %s under metID curation data of the tsv file', string(met(i)));
        end
        
        if ~ismember(alternativeID(i),'[]') && contains(alternativeID(i),'C')
            model.metNotes(idx_met(j)) = join([model.metNotes(idx_met(j)),'| alternative KEGG ID',alternativeID(i),'(PR #220)']);
        elseif ~ismember(alternativeID(i),'[]') && ~contains(alternativeID(i),'C')
            warning('Check for error in %s under metID curation data of the tsv file', string(met(i)));
        end
    end
end

%% Update metChEBIID

currentID = modifymetID(:,8);
newID = modifymetID(:,9);
alternativeID = modifymetID(:,10);

for i = 1:length(metNames)
    idx_met = find(ismember(modelR.metNames,metNames(i)));
    %Check for new ID, replace current ID/add new ID into model
    %if blank = remove current ID in model
    for j = 1:length(idx_met)
        if ~ismember(newID(i),'[]') && contains(newID(i),'CHEBI')
            model.metChEBIID(idx_met(j)) = newID(i);
            model.metNotes(idx_met(j)) = join([model.metNotes(idx_met(j)),'| ChEBI curated (PR #220)']);
        elseif ismember(newID(i),'Blank')
            model.metChEBIID(idx_met(j)) = {''};
            model.metNotes(idx_met(j)) = join([model.metNotes(idx_met(j)),'| ChEBI curated (PR #220)']);
        elseif ~ismember(newID(i),'[]') && ~contains(newID(i),'CHEBI')
            warning('Check for error in %s under metID curation data of the tsv file', string(met(i)));
        end
        
        if ~ismember(alternativeID(i),'[]') && contains(alternativeID(i),'CHEBI')
            model.metNotes(idx_met(j)) = join([model.metNotes(idx_met(j)),'| alternative ChEBI',alternativeID(i),'(PR #220)']);
        elseif ~ismember(alternativeID(i),'[]') && ~contains(alternativeID(i),'CHEBI')
            warning('Check for error in %s under metID curation data of the tsv file', string(met(i)));
        end
    end
end

%remove whitespace(s) when adding notes into rxnNotes
model.metNotes(:) = strtrim(model.metNotes(:));
% Remove leading '| ' in notes that were previously empty
model.rxnNotes = regexprep(model.rxnNotes,'^\| ','');
model.metNotes = regexprep(model.metNotes,'^\| ','');
%Save model
saveYeastModel(model);
cd modelCuration/
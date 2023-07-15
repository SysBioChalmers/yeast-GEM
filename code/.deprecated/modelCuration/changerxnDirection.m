% This script modifies reaction direction based on data from 
% rxnDirectionInfo.tsv
%
% rxnDirectionInfo.tsv is a compilation of deltaG, reversibility index and
% reaction direction information from various data sources, generated via
% checkrxnDirection.m. The tsv file is edited to add 3 new columns
% 'Action', 'New_model.lb' and 'Notes'
%
% Inputs: model and rxnDirectionInfo.tsv
%
% Cheng Wei Quan (Eiden), 2020-05-05

%Load model
cd ..
model = loadYeastModel;

%Load rxnDirectionInfo.tsv files
fid = fopen('../data/modelcuration/rxnDirectionInfo.tsv');
format = repmat('%s ',1,16);
format = strtrim(format);
temp = textscan(fid,format,'Delimiter','\t','HeaderLines',0);
for i = 1:length(temp)
    rxnDirectionInfo(:,i) = temp{i}; %use {} instead of () for cell array
end
commentLines = startsWith(rxnDirectionInfo(:,1),'#');
rxnDirectionInfo(commentLines,:) = [];
fclose(fid);

%Store required variables
rxns = rxnDirectionInfo(:,1);
[~,rxn_idx] = ismember(rxns,model.rxns);
action = rxnDirectionInfo(:,14);
new_lb = rxnDirectionInfo(:,15);

%Change rxnDirection/coefficient sign based on action
for i = 1:size(action,1)
    if rxn_idx(i) == 0
        warning('model.rxns %s cannot be found, check if it is an error', string(rxns(i)));
    elseif model.lb(rxn_idx(i)) == str2double(new_lb{i})
        %no change if model.lb already matches new_lb
    elseif strcmp(action(i),'Change to irreversible forward')
        model.lb(rxn_idx(i)) = str2double(new_lb{i});
        model.rxnNotes(rxn_idx(i)) = join([model.rxnNotes(rxn_idx(i)),'| rxnDirection curated (PR #227)']);
    elseif strcmp(action(i),'Change to irreversible forward and perform sign change for all metabolites in this reaction')
        model.lb(rxn_idx(i)) = str2double(new_lb{i});
        met_idx = find(model.S(:,rxn_idx(i)));
        model.S(met_idx,rxn_idx(i)) = -model.S(met_idx,rxn_idx(i));
        model.rxnNotes(rxn_idx(i)) = join([model.rxnNotes(rxn_idx(i)),'| rxnDirection curated (PR #227)']);
    end
end

%Remove whitespace(s) when adding notes
model.metNotes(:) = strtrim(model.metNotes(:));
model.rxnNotes(:) = strtrim(model.rxnNotes(:));

%Remove leading '| ' in notes that were previously empty
model.rxnNotes = regexprep(model.rxnNotes,'^\| ','');
model.metNotes = regexprep(model.metNotes,'^\| ','');

%Remove generic reaction r_4229, rxnName: Monocarboxylic acid amide amidohydrolase
[~,rxn_idx] = ismember('Monocarboxylic acid amide amidohydrolase',model.rxnNames);
if rxn_idx~=0
    model = removeRxns(model,model.rxns(rxn_idx)); 
end

%Save model
saveYeastModel(model);

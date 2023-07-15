% This script fixes unbalanced reaction in the model based on data from modMetsandSmatrix.tsv
% 
% modMetsandSmatrix.tsv includes details on changes to current metFormula
% and metCharges as well modifications of S matrix coefficient,
% after manual curation of selected unbalanced reactions in the model
%
% Inputs: model and modMetsandSmatrix.tsv
% Note: functions checkMassChargeBalance.m and changerxn.m are required
%

%Load model
cd ..
model = loadYeastModel;

%Load modMetsandSmatrix.tsv
fid = fopen('../data/modelCuration/modMetsandSmatrix.tsv');
format = repmat('%s ',1,14);
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
idx_mets = idx(1)+1:idx(2)-1;
idx_Smatrix = idx(2)+1:idx(3)-1;

%Separate data into various cell arrays
updatemets = curationfile(idx_mets,1:11);
updateSmatrix = curationfile(idx_Smatrix,1:9);

%Add 2 new metabolites to model 
metList = {'trans-4-hydroxy-L-proline [cytoplasm]';'2,3-dihydroxy-3-methylbutanoate [cytoplasm]'};
metFormula = {'C5H9NO3';'C5H9O4'};
metCharge = {0;-1};
metChEBIID = {'CHEBI:18072';'CHEBI:11424'};
metKEGGID = {'C01157';'C04039'};
metMetaNetXID = {'MNXM87584';'MNXM734'};

comps = split(metList, ' [');
comps = comps(:,2);
comps = strrep(comps,']','');
CONValldata = cat(2,model.compNames,model.comps);
[~,b] = ismember(comps,CONValldata(:,1));
comps = CONValldata(b,2);

for i = 1:length(metList)
    cd otherChanges/
    newID = getNewIndex(model.mets);
    cd ../
    mets(i) = strcat('s_',newID,'[',comps(i),']');
    model   = addMetabolite(model,char(mets(i)), ...
                            'metName',metList{i},'metFormula', ...
                            metFormula{i},'Charge',metCharge{i}, ...
                            'ChEBIID',metChEBIID{i},'KEGGId',metKEGGID{i});
    %Manually add MetaNetXID as addMetabolite does not include it
    met_idx = find(ismember(model.mets,mets(i)));
    model.metMetaNetXID(met_idx) = metMetaNetXID(i);
end

cd modelCuration/

%Modify the following rxns:
%r_0687: L-proline is replaced by trans-4-hydroxy-L-proline
model = changerxn(model, 'r_0687', '1-pyrroline-3-hydroxy-5-carboxylic&acid [cytoplasm] + 2 H+ [cytoplasm] + NADPH [cytoplasm]  -> trans-4-hydroxy-L-proline [cytoplasm] + NADP(+) [cytoplasm]');
%r_4577: 3-hydroxy-3-methyl-2-oxobutanoate is replaced by 2,3-dihydroxy-3-methylbutanoate
model = changerxn(model, 'r_4577', '3-methyl-2-oxobutanoate [cytoplasm] + H2O [cytoplasm]  -> 2,3-dihydroxy-3-methylbutanoate [cytoplasm]');
cd ../missingFields/
model = addSBOterms(model); %Add SBO terms
cd ../modelCuration/
model = rmfield(model,'grRules'); %remove field 'grRules' in model

%Correction of metFormula and metCharges to balance equation
met = updatemets(:,1);
[~,idx_met] = ismember(met,model.mets);
modelR = ravenCobraWrapper(model);
metNames = modelR.metNames(idx_met);
currentFormula = updatemets(:,2);
newFormula = updatemets(:,3);
currentCharge = str2double(updatemets(:,4));
newCharge = str2double(updatemets(:,5));
MNXNotes = updatemets(:,11);
metResults{size(updatemets,1)*5,6} = [];
metResults(1,:) = [{'Metabolites'},{'Unbalanced reaction(s) before modification (UB)'},...
                   {'Unbalanced reaction after modification (UA)'},{'Number of UB'},...
                   {'Number of UA'},{'Number of balanced reactions'}];

for i = 1:length(metNames)
    idx_met = find(ismember(modelR.metNames,metNames(i)));
    currentmetFormula = model.metFormulas(idx_met);
    currentmetCharges = model.metCharges(idx_met);
    
    for j = 1:length(idx_met)
        %Preliminary check of metBalance/metCharges before change
        idx_rxn = find(model.S(idx_met(j),:));
        rxn = model.rxns(idx_rxn);
        MassChargeresults = checkMassChargeBalance(model,rxn);
            
        if isequal(currentmetFormula(j),newFormula(i))
            %no action required
        elseif ~isequal(currentmetFormula(j),currentFormula(i)) && ismember(currentFormula(i),'[]') || isequal(currentmetFormula(j),currentFormula(i)) && ~ismember(newFormula(i),'[]')
            model.metFormulas(idx_met(j)) = newFormula(i);
            if ~contains(model.metNotes(idx_met(j)),'| metFormula curated (PR #222)')
                model.metNotes(idx_met(j)) = join([model.metNotes(idx_met(j)),'| metFormula curated (PR #222)']);
            end
        elseif ~isequal(currentmetFormula(j),currentFormula(i))
            warning('error with metFormula matching in modMetsandSmatrix.tsv for %s, idx_met: %d', string(metNames(i)), idx_met(j));
        end
        
        if isequal(currentmetCharges(j),newCharge(i))
            %no action required
        elseif (isequal(currentmetCharges(j),currentCharge(i)) && ~isnan(newCharge(i))) || (ismissing(currentmetCharges(j)) && isnan(currentCharge(i)) && ~isnan(newCharge(i)))
            model.metCharges(idx_met(j)) = newCharge(i);
            if ~contains(model.metNotes(idx_met(j)),'| metCharge curated (PR #222)')
                model.metNotes(idx_met(j)) = join([model.metNotes(idx_met(j)),'| metCharge curated (PR #222)']);
            end
        elseif ~isequal(currentmetCharges(j),currentCharge(i)) && ~ismissing(currentmetCharges(j))
            warning('error with metCharges matching in modMetsandSmatrix.tsv for %s, idx_met: %d', string(metNames(i)), idx_met(j));
        end
        
        %Find rxns which are unbalanced before change
        temp_rxn = MassChargeresults(:,1);
        temp_notes = MassChargeresults(:,2);
        temp_idx = find(~ismember(temp_notes,'pass'));
        unbalancedrxnBefore = temp_rxn(temp_idx);
        unbalancedNotes = temp_notes(temp_idx);
        unbalancedBefore = length(find(~ismember(temp_notes,'pass')));
        
        arrayidx = find(cellfun('isempty', metResults),1);
        metResults(arrayidx,1) = cellstr(metNames(i));
        metResults(arrayidx,2) = {[unbalancedrxnBefore,unbalancedNotes]};
        metResults(arrayidx,4) = num2cell(unbalancedBefore);
    end
end

%Modification of H+/H2O/other metabolite coefficient to balance equation
met2 = updateSmatrix(:,1);
rxn2 = updateSmatrix(:,2);
[~,idx_met2] = ismember(met2,model.mets);
[~,idx_rxn2] = ismember(rxn2,model.rxns);
currentCoef = str2double(updateSmatrix(:,3));
newCoef = str2double(updateSmatrix(:,4));

for i = 1:length(updateSmatrix)
    currentmetCoef = model.S(idx_met2(i),idx_rxn2(i));
    if isequal(currentmetCoef,currentCoef(i)) && ~isnan(newCoef(i))
        model.S(idx_met2(i),idx_rxn2(i)) = newCoef(i);
        model.rxnNotes(idx_rxn2(i)) = join([model.rxnNotes(idx_rxn2(i)),'| model.S(',...
                cellstr(string(idx_met2(i))),',',cellstr(string(idx_rxn2(i))),') curated (PR #222)']);
        model.rxnNotes(idx_rxn2(i)) = strrep(model.rxnNotes(idx_rxn2(i)),'( ','(');
        model.rxnNotes(idx_rxn2(i)) = strrep(model.rxnNotes(idx_rxn2(i)),' )',')');
        model.rxnNotes(idx_rxn2(i)) = strrep(model.rxnNotes(idx_rxn2(i)),' , ',',');
    elseif ~isequal(currentmetCoef,currentCoef(i))
        warning('error with metCoef matching in modMetsandSmatrix.tsv for %s, idx_rxn: %d', string(rxn2(i)), idx_rxn2(i));
    end
end 

%remove whitespace(s) when adding notes
model.metNotes(:) = strtrim(model.metNotes(:));
model.rxnNotes(:) = strtrim(model.rxnNotes(:));

%Check of metFormula/metCharges after change (for updatemets)
for i = 1:length(metNames)
    idx_met = find(ismember(modelR.metNames,metNames(i)));
    for j = 1:length(idx_met)
        idx_rxn = find(model.S(idx_met(j),:));
        rxn = model.rxns(idx_rxn);
        MassChargeresults2 = checkMassChargeBalance(model,rxn);

        %Find rxns which are unbalanced after change
        temp_rxn2 = MassChargeresults2(:,1);
        temp_notes2 = MassChargeresults2(:,2);
        temp_idx2 = find(~ismember(temp_notes2,'pass'));
        unbalancedrxnAfter = temp_rxn2(temp_idx2);
        unbalancedNotes2 = temp_notes2(temp_idx2);
        unbalancedAfter = length(find(~ismember(temp_notes2,'pass')));
        
        if isempty(unbalancedrxnAfter)
            arrayidx = find(cellfun('isempty', metResults(:,3)),1);
            metResults(arrayidx,3) = {'NIL'};
            metResults(arrayidx,5) = {0};
            netChange = str2double(string((metResults(arrayidx,4)))) - unbalancedAfter;
            metResults(arrayidx,6) = num2cell(netChange);
        else
            arrayidx = find(cellfun('isempty', metResults(:,3)),1);
            metResults(arrayidx,3) = {[unbalancedrxnAfter,unbalancedNotes2]};
            metResults(arrayidx,5) = num2cell(unbalancedAfter);
            netChange = str2double(string((metResults(arrayidx,4)))) - unbalancedAfter;
            metResults(arrayidx,6) = num2cell(netChange);
        end
    end
end

%Check of metBalance/metCharges after change (for updateSmatrix)
metResults2{size(updateSmatrix,1),5} = [];
for i = 1:size(updateSmatrix,1)
    rxn = updateSmatrix(i,2);
    MassChargeresults2 = checkMassChargeBalance(model,rxn);
    metResults2(i,:) = MassChargeresults2;
end

%Remove leading '| ' in notes that were previously empty
model.rxnNotes = regexprep(model.rxnNotes,'^\| ','');
model.metNotes = regexprep(model.metNotes,'^\| ','');

%Save model
cd ..
saveYeastModel(model);

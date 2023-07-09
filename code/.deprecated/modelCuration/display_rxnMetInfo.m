% This function is to output relevant metabolite data of a single reaction
% to facilitate manual curation
%
% Input: model and reaction identifier in model.rxns (e.g. 'r_0001')
% Output: editFormula (n x 9 cell array, n = number of mets in rxn)
%         
%         All mets info are displayed in the following order:
%         model.mets, (met_idx, rxn_idx), model.metNames, 
%         model.metMetaNetXID,model.metKEGGID, model.metChEBIID,
%         model.S, model.metCharges
%
% NOTE: reac_prop.tsv and chem_prop.tsv will be downloaded in current 
%       directory when using this function for the first time, if they
%       cannot be not found - will increase the run time
%
% Usage: editFormula = display_rxnMetInfo(model,rxn_num)
%

function editFormula = display_rxnMetInfo(model,rxn_num)
%Check for reac_prop.tsv in current directory
%If not available, file will be downloaded
downloadMNXdb('reac_prop',pwd)

%Load reac_prop.tsv files from MNX database
fid2 = fopen('reac_prop.tsv');
format = repmat('%s ',1,6);
format = strtrim(format);
rxn_temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',0);
for i = 1:length(rxn_temp)
    MNXreacprop(:,i) = rxn_temp{i};
end
commentLines = startsWith(MNXreacprop(:,1),'#');
MNXreacprop(commentLines,:) = []; % MNX_ID Equation Description Balance EC Source
fclose(fid2);

%Check for chem_prop.tsv in current directory
%If not available, file will be downloaded
downloadMNXdb('chem_prop',pwd)

%Load chem_prop.tsv files from MNX database
fid2 = fopen('chem_prop.tsv');
format = repmat('%s ',1,9);
format = strtrim(format);
met_temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',0);
for i = 1:length(met_temp)
    MNXchemprop(:,i) = met_temp{i};
end
commentLines = startsWith(MNXchemprop(:,1),'#');
MNXchemprop(commentLines,:) = []; % MNX_ID Description Formula Charge Mass InChI SMILES Source InChIKey
fclose(fid2);

%find the formula which causes equation to be unbalanced
bs  = getElementalBalance(model);
idx = find(bs.balanceStatus==0);
out = model.rxns(idx);
rxn_idx = find(ismember(model.rxns,rxn_num));

%find element difference
dif = bs.leftComp(rxn_idx,:) - bs.rightComp(rxn_idx,:);
mets = find(~dif==0);
difference = '';
for j=1:length(mets)
    difference = strcat(difference,bs.elements.abbrevs(mets(j)));
    difference = strcat(difference,num2str(dif(mets(j))));
end

%Find metFormulas in model
met_idx = find(model.S(:,rxn_idx));
metNames = model.metNames(met_idx);
mets = model.mets(met_idx);
metFormulas = model.metFormulas(met_idx);
MNXmetID = model.metMetaNetXID(met_idx);
KEGGmetID = model.metKEGGID(met_idx);
ChEBImetID = model.metChEBIID(met_idx);
SmatrixCoef{length(met_idx),1} = [];
for i = 1:length(met_idx)
    SmatrixCoef(i) = num2cell(model.S(met_idx(i),rxn_idx));
end
metCharges = model.metCharges(met_idx);
arrayrows = max([length(MNXmetID),length(metFormulas),length(metCharges)]);

editFormula{arrayrows,7} = [];
editFormula(1:length(mets),1) = mets;
for i = 1:length(met_idx)
    editFormula(i,2) = {sprintf('%s, %s',string(met_idx(i)),string(rxn_idx))};
end
editFormula(1:length(MNXmetID),3) = metNames;
editFormula(1:length(MNXmetID),4) = MNXmetID;
editFormula(1:length(KEGGmetID),5) = KEGGmetID;
editFormula(1:length(ChEBImetID),6) = ChEBImetID;
editFormula(1:length(metFormulas),7) = metFormulas;
editFormula(1:length(SmatrixCoef),8) = SmatrixCoef;
editFormula(1:length(metCharges),9) = num2cell(metCharges);

editFormula(cellfun('isempty',editFormula)) = {0};
editFormula = string(editFormula);

end

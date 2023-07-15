% This script maps metabolite IDs of different databases to find new MNXMID(s)
%
% Inputs: model
%
% New MNXMID(s) will only be added if 1 of the following 3 criteria is met:
%   a) Both metKEGGID and metChEBIID is mapped to the same MNXMID
%   b) Only metKEGGID mapped to MNXMID, metFormula of new MNXMID matches model.metFormula 
%   c) Only metChEBIID mapped to MNXMID, metFormula of new MNXMID matches model.metFormula 
%
% Inputs: model
%

%load model
cd ..
model = loadYeastModel;

%Check for chem_prop.tsv in current directory
%If not available, file will be downloaded
downloadMNXdb('chem_prop',pwd)

%Load MNXchem_prop.tsv file containing data on all compounds from MetaNetX database
fid = fopen('chem_prop.tsv');
format = repmat('%s ',1,9);
format = strtrim(format);
met_temp = textscan(fid,format,'Delimiter','\t','HeaderLines',0);
for i = 1:length(met_temp)
    MNXchem_prop(:,i) = met_temp{i};
end
commentLines = startsWith(MNXchem_prop(:,1),'#');
MNXchem_prop(commentLines,:) = []; %MNX_ID Description Formula Charge Mass InChI SMILES Source InChIKey
fclose(fid);

%map metCheBIID to metMetaNetXID via function mapIDsViaMNXref
cd missingFields
temp = replace(model.metChEBIID,'CHEBI:','');
xref_metMetaNetX = mapIDsViaMNXref('mets',temp,'ChEBI','MetaNetX');
xref_metMetaNetX(:,2) = model.mets; %match MNXMID with model.mets
xref_metMetaNetX(:,3) = model.metChEBIID; %match MNXMID with model.metChEBIID
xref_metMetaNetX(:,4) = model.metMetaNetXID; %match MNXMID with model.metMetaNetXID

%only retain mapped MNXMID if model.metMetaNetXID is empty
empties = find(cellfun('isempty',xref_metMetaNetX(:,1)));
xref_metMetaNetX(empties,:) = [];
empties = find(~cellfun('isempty',xref_metMetaNetX(:,4)));
xref_metMetaNetX(empties,:) = [];

%map metKEGGID to metMetaNetXID via function mapIDsViaMNXref
xref_metMetaNetX_2 = mapIDsViaMNXref('mets',model.metKEGGID,'KEGG','MetaNetX');
xref_metMetaNetX_2(:,2) = model.mets; %match MNXMID with model.mets
xref_metMetaNetX_2(:,3) = model.metKEGGID; %match MNXMID with model.metKEGGID
xref_metMetaNetX_2(:,4) = model.metMetaNetXID; %match MNXMID with model.metMetaNetXID

%only retain mapped MNXMID if model.metMetaNetXID is empty
empties = find(cellfun('isempty',xref_metMetaNetX_2(:,1)));
xref_metMetaNetX_2(empties,:) = [];
empties = find(~cellfun('isempty',xref_metMetaNetX_2(:,4)));
xref_metMetaNetX_2(empties,:) = [];

%check for cases in which MNXMID mapped via both metChEBIID and metKEGGID
%compile data required via indexing
[~,idx] = ismember(xref_metMetaNetX_2(:,2),xref_metMetaNetX(:,2));
idx = idx(idx~=0);
matchKEGGChEBI(:,1) = xref_metMetaNetX(idx,2); %model.mets
matchKEGGChEBI(:,3) = xref_metMetaNetX(idx,1); %MNXMID mapped via metChEBIID
matchKEGGChEBI(:,4) = xref_metMetaNetX(idx,3); %model.metChEBIID

[~,idx2] = ismember(matchKEGGChEBI(:,1),model.mets);
idx2 = idx2(idx2~=0);
matchKEGGChEBI(:,2) = model.metFormulas(idx2); %model.metFormulas

[~,idx3] = ismember(matchKEGGChEBI(:,3),MNXchem_prop(:,1));
matchKEGGChEBI(:,5) = {''};
matchKEGGChEBI(idx3~=0,5) = MNXchem_prop(idx3(idx3~=0),3); %metFormula of MNXMID mapped via metChEBIID

[~,idx4] = ismember(xref_metMetaNetX(:,2),xref_metMetaNetX_2(:,2));
idx4 = idx4(idx4~=0);
matchKEGGChEBI(:,6) = xref_metMetaNetX_2(idx4,1); %MNXMID mapped via metKEGGID
matchKEGGChEBI(:,7) = xref_metMetaNetX_2(idx4,3); %model.metKEGGID

[~,idx5] = ismember(matchKEGGChEBI(:,6),MNXchem_prop(:,1));
idx5 = idx5(idx5~=0);
matchKEGGChEBI(:,8) = MNXchem_prop(idx5,3); %metFormula of MNXMID mapped via metKEGGID

for i = 1:size(matchKEGGChEBI,1)
    if ismember(matchKEGGChEBI(i,3),matchKEGGChEBI(i,6))
        %add MNXMID into model if both metKEGGID and metChEBIID is mapped to the same MNXMID
        model.metMetaNetXID(idx2(i)) = matchKEGGChEBI(i,3);
        matchKEGGChEBI(i,9) = join(['added',model.metMetaNetXID(idx2(i)),'into model']);
    else
        if ismember(matchKEGGChEBI(i,2),matchKEGGChEBI(i,5))
            %add MNXMID into model if MNXMID mapped via metChEBIID has same metFormula as model.metFormula
            model.metMetaNetXID(idx2(i)) = matchKEGGChEBI(i,3);
            matchKEGGChEBI(i,9) = join(['added',model.metMetaNetXID(idx2(i)),'into model']);
        elseif ismember(matchKEGGChEBI(i,2),matchKEGGChEBI(i,8))
            %add MNXMID into model if MNXMID mapped via metKEGGID has same metFormula as model.metFormula
            model.metMetaNetXID(idx2(i)) = matchKEGGChEBI(i,3);
            matchKEGGChEBI(i,9) = join(['added',model.metMetaNetXID(idx2(i)),'into model']);
        else
            matchKEGGChEBI(i,9) = {'MNXMID not added'};
        end
    end      
end

%check for cases in which MNXMID mapped via metChEBIID only
%compile data required via indexing
idx = find(~ismember(xref_metMetaNetX(:,2),xref_metMetaNetX_2(:,2)));
matchChEBI(:,1) = xref_metMetaNetX(idx,2); %model.mets
matchChEBI(:,3) = xref_metMetaNetX(idx,1); %MNXMID mapped via metChEBIID
matchChEBI(:,4) = xref_metMetaNetX(idx,3); %model.metChEBIID

[~,idx2] = ismember(matchChEBI(:,1),model.mets);
idx2 = idx2(idx2~=0);
matchChEBI(:,2) = model.metFormulas(idx2); %model.metFormulas

[~,idx3] = ismember(matchChEBI(:,3),MNXchem_prop(:,1));
matchChEBI(:,5) = {''};
matchChEBI(idx3~=0,5) = MNXchem_prop(idx3(idx3~=0),3); %metFormula of MNXMID mapped via metChEBIID

for i = 1:size(matchChEBI,1)
    if ismember(matchChEBI(i,2),matchChEBI(i,5))
        %add MNXMID into model if MNXMID mapped via metChEBIID has same metFormula as model.metFormula
        model.metMetaNetXID(idx2(i)) = matchChEBI(i,3);
        matchChEBI(i,6) = join(['added',model.metMetaNetXID(idx2(i)),'into model']);
    else
        matchChEBI(i,6) = {'MNXMID not added'};
    end
end

%check for cases in which MNXMID mapped via metKEGGID only
%compile data required via indexing
idx = find(~ismember(xref_metMetaNetX_2(:,2),xref_metMetaNetX(:,2)));
matchKEGG(:,1) = xref_metMetaNetX_2(idx,2); %model.mets
matchKEGG(:,3) = xref_metMetaNetX_2(idx,1); %MNXMID mapped via metKEGGID
matchKEGG(:,4) = xref_metMetaNetX_2(idx,3); %model.metKEGGID

[~,idx2] = ismember(matchKEGG(:,1),model.mets);
idx2 = idx2(idx2~=0);
matchKEGG(:,2) = model.metFormulas(idx2); %model.metFormulas

[~,idx3] = ismember(matchKEGG(:,3),MNXchem_prop(:,1));
matchKEGG(:,5) = {''};
matchKEGG(idx3~=0,5) = MNXchem_prop(idx3(idx3~=0),3); %metFormula of MNXMID mapped via metKEGGID

for i = 1:size(matchKEGG,1)
    if ismember(matchKEGG(i,2),matchKEGG(i,5))
        %add MNXMID into model if MNXMID mapped via metChEBIID has same metFormula as model.metFormula
        model.metMetaNetXID(idx2(i)) = matchKEGG(i,3);
        matchKEGG(i,6) = join(['added',model.metMetaNetXID(idx2(i)),'into model']);
    else
        matchKEGG(i,6) = {'MNXMID not added'};
    end
end

%Save model
cd ..
fileDir = dir('chem_prop.tsv');
if ~isempty(fileDir)
    delete('chem_prop.tsv');
end
saveYeastModel(model);
cd modelCuration/
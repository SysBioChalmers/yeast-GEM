% This script maps reaction IDs and metabolite IDs of different databases
% to find new rxnKEGGID(s) and metKEGGID(s)
%
% Inputs: model
%
% getRxnsFromMetaCyc and getMetsFromMetaCyc are functions from RAVEN
%
% reactions.tsv and compounds.tsv from modelSEED database is downloaded directly from GitHub repository ModelSEEDDatabase
% link: github.com/ModelSEED/ModelSEEDDatabase
%

%load model
cd ..
model = loadYeastModel;

%% rxnKEGGID
%map rxnMetaNetXID to rxnKEGGID via function mapIDsViaMNXref
cd missingFields
xref_rxnKEGG = mapIDsViaMNXref('rxns',model.rxnMetaNetXID,'MetaNetX','KEGG');
xref_rxnKEGG(:,2) = model.rxns; %match rxnKEGGID with model.rxns
xref_rxnKEGG(:,3) = model.rxnMetaNetXID; %match rxnKEGGID with model.rxnMetaNetXID
xref_rxnKEGG(:,4) = model.rxnKEGGID; %match rxnKEGGID with model.rxnKEGGID

%only retain mapped rxnKEGGID if model.rxnKEGGID is empty
empties = find(cellfun('isempty',xref_rxnKEGG(:,1)));
xref_rxnKEGG(empties,:) = [];
empties = find(~cellfun('isempty',xref_rxnKEGG(:,4)));
xref_rxnKEGG(empties,:) = [];

%add rxnKEGGID into model
[~,idx] = ismember(xref_rxnKEGG(:,2),model.rxns);
idx = idx(idx~=0);

for i = 1:length(idx)
    model.rxnKEGGID(idx(i)) = xref_rxnKEGG(i,1);
end

%correction for MNXref error in MNXR96123/R00124#2
idx = find(ismember(model.rxnMetaNetXID,'MNXR96123'));
model.rxnKEGGID(idx) = {'R00124'};

%map rxnMetaNetXID to rxnMetaCycID via function mapIDsViaMNXref, then map rxnMetaCycID to rxnKEGGID
xref_rxnMetaCyc = mapIDsViaMNXref('rxns',model.rxnMetaNetXID,'MetaNetX','MetaCyc');
xref_rxnMetaCyc(:,2) = model.rxns; %match rxnMetaCycID with model.rxns
xref_rxnMetaCyc(:,3) = model.rxnMetaNetXID; %match rxnMetaCycID with model.rxnMetaNetXID
xref_rxnMetaCyc(:,4) = model.rxnKEGGID; %match rxnMetaCycID with model.rxnKEGGID
empties = find(cellfun('isempty',xref_rxnMetaCyc(:,1)));
xref_rxnMetaCyc(empties,:) = [];

%load metaCycRxns.mat via getRxnsFromMetaCyc from RAVEN directory
MetaCyc_rxnInfo = [];
try
    MetaCyc_rxnInfo = getRxnsFromMetaCyc;
catch
    warning('RAVEN repository is not cloned or added to MATLAB path, unable to use function getRxnsFromMetaCyc');
end

%Input rxnKEGGID into xref array
if ~isempty(MetaCyc_rxnInfo) %check if MetaCyc_rxnInfo is generated successfully
    for i = 1:length(xref_rxnMetaCyc)
        if isempty(xref_rxnMetaCyc{i,4})
            idx = find(ismember(MetaCyc_rxnInfo.rxns,xref_rxnMetaCyc(i,1)));
            if idx~=0
                if ~isempty(MetaCyc_rxnInfo.rxnMiriams{idx})
                    temp = cell2mat(MetaCyc_rxnInfo.rxnMiriams(idx));
                    idx2 = find(ismember(temp.name,'kegg.reaction'));
                    if ~isempty(idx2)
                        xref_rxnMetaCyc(i,5) = temp.value(idx2);
                    end
                end
            end
        end
    end
    
    if size(xref_rxnMetaCyc,2) == 5
        %only retain mapped rxnKEGGID if model.rxnKEGGID is empty
        empties = find(cellfun('isempty',xref_rxnMetaCyc(:,5)));
        xref_rxnMetaCyc(empties,:) = [];
        empties = find(~cellfun('isempty',xref_rxnMetaCyc(:,4)));
        xref_rxnMetaCyc(empties,:) = [];
        
        %add rxnKEGGID into model
        [~,idx] = ismember(xref_rxnMetaCyc(:,2),model.rxns);
        idx = idx(idx~=0);
        
        for i = 1:length(idx)
            model.rxnKEGGID(idx(i)) = xref_rxnMetaCyc(i,5);
        end
    else
        xref_rxnMetaCyc(:) = []; %erase all cell array, since no KEGGID matched and added
    end
else
    warning('ID mapping with MetaCyc IDs unsuccessful as function getMetsFromMetaCyc cannot be used to generate MetaCyc_metInfo');
end

%map rxnMetaNetXID to rxnSEEDID via function mapIDsViaMNXref, then map rxnSEEDID to rxnKEGGID
xref_rxnSEED = mapIDsViaMNXref('rxns',model.rxnMetaNetXID,'MetaNetX','SEED');
xref_rxnSEED(:,2) = model.rxns; %match rxnSEEDID with model.rxns
xref_rxnSEED(:,3) = model.rxnMetaNetXID; %match rxnSEEDID with model.rxnMetaNetXID
xref_rxnSEED(:,4) = model.rxnKEGGID; %match rxnSEEDID with model.rxnKEGGID
empties = find(cellfun('isempty',xref_rxnSEED(:,1)));
xref_rxnSEED(empties,:) = [];

%download reactions.tsv from GitHub repository ModelSEED/ModelSEEDDatabase
try
    websave('reactions.tsv','https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/dev/Biochemistry/reactions.tsv');
catch
    warning('reactions.tsv was not successfully downloaded, check if directory for reactions.tsv has changed on github.com/ModelSEED/ModelSEEDDatabase/Biochemistry');
end
fileDir = dir('reactions.tsv');

if ~isempty(fileDir) %check if reactions.tsv is present in current directory
    %load reactions.tsv
    fid2 = fopen('reactions.tsv');
    format = repmat('%s ',1,22);
    format = strtrim(format);
    rxn_temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',0);
    for i = 1:length(rxn_temp)
        SEED_rxnInfo(:,i) = rxn_temp{i};
    end
    commentLines = startsWith(SEED_rxnInfo(:,1),'#');
    SEED_rxnInfo(commentLines,:) = [];
    fclose(fid2);
    
    %Input rxnKEGGID into xref array
    for i = 1:length(xref_rxnSEED)
        if isempty(xref_rxnSEED{i,4})
            idx = find(ismember(SEED_rxnInfo(:,1),xref_rxnSEED(i,1)));
            if idx~=0
                if ~isempty(SEED_rxnInfo{idx,13}) && contains(SEED_rxnInfo(idx,13),'KEGG')
                    temp_ID = extractBetween(SEED_rxnInfo(idx,13),'KEGG: ','|');
                    xref_rxnSEED(i,5) = temp_ID;
                end
            end
        end
    end
    
    if size(xref_rxnSEED,2) == 5
        %only retain mapped rxnKEGGID if model.rxnKEGGID is empty
        empties = find(cellfun('isempty',xref_rxnSEED(:,5)));
        xref_rxnSEED(empties,:) = [];
        empties = find(~cellfun('isempty',xref_rxnSEED(:,4)));
        xref_rxnSEED(empties,:) = [];
        
        %add rxnKEGGID into model
        [~,idx] = ismember(xref_rxnSEED(:,2),model.rxns);
        idx = idx(idx~=0);
        
        for i = 1:length(idx)
            model.rxnKEGGID(idx(i)) = xref_rxnSEED(i,5);
        end
    else
        xref_rxnSEED(:) = []; %erase all cell array, since no KEGGID matched and added
    end
    delete('reactions.tsv');
else
    warning('ID mapping with modelSEED IDs unsuccessful as reactions.tsv is not found');
end

%% metKEGGID
%map metMetaNetXID to metKEGGID via function mapIDsViaMNXref
xref_metKEGG = mapIDsViaMNXref('mets',model.metMetaNetXID,'MetaNetX','KEGG');
xref_metKEGG(:,2) = model.mets; %match metKEGGID with model.mets
xref_metKEGG(:,3) = model.metMetaNetXID; %match metKEGGID with model.metMetaNetXID
xref_metKEGG(:,4) = model.metKEGGID; %match metKEGGID with model.metKEGGID

%only retain mapped metKEGGID if model.metKEGGID is empty
empties = find(cellfun('isempty',xref_metKEGG(:,1)));
xref_metKEGG(empties,:) = [];
empties = find(~cellfun('isempty',xref_metKEGG(:,4)));
xref_metKEGG(empties,:) = [];

%add metKEGGID into model
[~,idx] = ismember(xref_metKEGG(:,2),model.mets);
idx = idx(idx~=0);

for i = 1:length(idx)
    if ~contains(xref_metKEGG(i,1),'G') %metKEGGID containing G not added, since it is not compatible with current SBML format
        model.metKEGGID(idx(i)) = xref_metKEGG(i,1);
    end
end

%map metMetaNetXID to metMetaCycID via function mapIDsViaMNXref, then map metMetaCycID to metKEGGID
xref_metMetaCyc = mapIDsViaMNXref('mets',model.metMetaNetXID,'MetaNetX','MetaCyc');
xref_metMetaCyc(:,2) = model.mets; %match metMetaCycID with model.mets
xref_metMetaCyc(:,3) = model.metMetaNetXID; %match metMetaCycID with model.metMetaNetXID
xref_metMetaCyc(:,4) = model.metKEGGID; %match metMetaCycID with model.metKEGGID
empties = find(cellfun('isempty',xref_metMetaCyc(:,1)));
xref_metMetaCyc(empties,:) = [];

%load metaCycMets.mat via getMetsFromMetaCyc from RAVEN
MetaCyc_metInfo = [];
try
    MetaCyc_metInfo = getMetsFromMetaCyc;
catch
    warning('RAVEN repository is not cloned or added to MATLAB path, unable to use function getMetsFromMetaCyc');
end

%Input metKEGGID into xref array
if ~isempty(MetaCyc_metInfo) %check if MetaCyc_metInfo is generated successfully
    for i = 1:length(xref_metMetaCyc)
        if isempty(xref_metMetaCyc{i,4})
            idx = find(ismember(MetaCyc_metInfo.mets,xref_metMetaCyc(i,1)));
            if idx~=0
                if ~isempty(MetaCyc_metInfo.keggid{idx})
                    xref_metMetaCyc(i,5) = MetaCyc_metInfo.keggid(idx);
                end
            end
        end
    end
    
    if size(xref_metMetaCyc,2) == 5
        %only retain mapped metKEGGID if model.metKEGGID is empty
        empties = find(cellfun('isempty',xref_metMetaCyc(:,5)));
        xref_metMetaCyc(empties,:) = [];
        empties = find(~cellfun('isempty',xref_metMetaCyc(:,4)));
        xref_metMetaCyc(empties,:) = [];
        
        %add metKEGGID into model
        [~,idx] = ismember(xref_metMetaCyc(:,2),model.mets);
        idx = idx(idx~=0);
        
        for i = 1:length(idx)
            if ~contains(xref_metMetaCyc(i,5),'G') %metKEGGID containing G not added, since it is not compatible with current SBML format
                model.metKEGGID(idx(i)) = xref_metMetaCyc(i,5);
            end
        end
    else
        xref_metMetaCyc(:) = []; %erase all cell array, since no KEGGID matched and added
    end
else
    warning('ID mapping with MetaCyc IDs unsuccessful as function getMetsFromMetaCyc cannot be used to generate MetaCyc_metInfo');
end

%map metMetaNetXID to metSEEDID via function mapIDsViaMNXref, then map metSEEDID to metKEGGID
xref_metSEED = mapIDsViaMNXref('mets',model.metMetaNetXID,'MetaNetX','SEED');
xref_metSEED(:,2) = model.mets; %match metSEEDID with model.mets
xref_metSEED(:,3) = model.metMetaNetXID; %match metSEEDID with model.metMetaNetXID
xref_metSEED(:,4) = model.metKEGGID; %match metSEEDID with model.metKEGGID
empties = find(cellfun('isempty',xref_metSEED(:,1)));
xref_metSEED(empties,:) = [];

%download reactions.tsv from GitHub repository ModelSEED/ModelSEEDDatabase
try
    websave('compounds.tsv','https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/dev/Biochemistry/compounds.tsv');
catch
    warning('compounds.tsv was not successfully downloaded, check if directory for compounds.tsv has changed on github.com/ModelSEED/ModelSEEDDatabase/Biochemistry');
end
fileDir = dir('compounds.tsv');

if ~isempty(fileDir) %check if compounds.tsv is present in current directory
    %load compounds.tsv
    fid2 = fopen('compounds.tsv');
    format = repmat('%s ',1,21);
    format = strtrim(format);
    met_temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',0);
    for i = 1:length(met_temp)
        SEED_metInfo(:,i) = met_temp{i};
    end
    commentLines = startsWith(SEED_metInfo(:,1),'#');
    SEED_metInfo(commentLines,:) = [];
    fclose(fid2);
    
    %Input metKEGGID into xref array
    for i = 1:length(xref_metSEED)
        if isempty(xref_metSEED{i,4})
            idx = find(ismember(SEED_metInfo(:,1),xref_metSEED(i,1)));
            if idx~=0
                if ~isempty(SEED_metInfo{idx,19}) && contains(SEED_metInfo(idx,19),'KEGG')
                    temp_ID = extractBetween(SEED_metInfo(idx,19),'KEGG: ','|');
                    if isempty(temp_ID) %KEGGID may be present at the end
                        temp_ID = extractAfter(SEED_metInfo(idx,19),'KEGG: ');
                        temp_ID = erase(temp_ID,'"');
                    end
                    if contains(temp_ID,';') %if more than 1 ID recorded, take the 1st one by default
                        temp_ID = split(temp_ID,';',2);
                        temp_ID = temp_ID(1,1);
                        warning('check %s, may contain more than 1 metKEGGID',string(SEED_metInfo(idx,1)));
                    end
                    xref_metSEED(i,5) = temp_ID;
                end
            end
        end
    end
    
    if size(xref_metSEED,2) == 5
        %only retain mapped metKEGGID if model.metKEGGID is empty
        empties = find(cellfun('isempty',xref_metSEED(:,5)));
        xref_metSEED(empties,:) = [];
        empties = find(~cellfun('isempty',xref_metSEED(:,4)));
        xref_metSEED(empties,:) = [];
        
        %add metKEGGID into model
        [~,idx] = ismember(xref_metSEED(:,2),model.mets);
        idx = idx(idx~=0);
        
        for i = 1:length(idx)
            if ~contains(xref_metSEED(i,5),'G') %metKEGGID containing G not added, since it is not compatible with current SBML format
                model.metKEGGID(idx(i)) = xref_metSEED(i,5);
            end
        end
    else
        xref_metSEED(:) = []; %erase all cell array, since no KEGGID matched and added
    end
    delete('compounds.tsv');
else
    warning('ID mapping with modelSEED IDs unsuccessful as compounds.tsv is not found');
end

%Save model
cd ..
saveYeastModel(model);
cd modelCuration/
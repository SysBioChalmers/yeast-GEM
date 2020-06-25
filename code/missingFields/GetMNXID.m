function model = GetMNXID(model,type,level)
% GetMNXID
%   mapping the model into metanetxID
%   In this function, we call the function 'mapIDsViaMNXref.m' from /ComplementaryScripts/missingFields, which is orginiated from
%   HMR3, and the part 'Load MNXref data structure' is modified to read the reference data from RAVEN.
%
%   model                   a model structure;please load yeast model using function 'loadYeastModel.m' from
%                           yeastGEM github repository.
%   type                    type = 'mets' or 'rxns'
%   level                   level = 1:only compares MNXid FROM kegg and CHEBI
%                           level =2 compares MNX from yeast7.6MNX model and CHEBI and KEGG
%                           (default =1)
%
%   Usage: model = GetMNXID(model,type,level)
%

if nargin<3
    level = 1;
end

if strcmpi(type,'mets')
    %mapping KEGGID to metanet ID
    query1 = model.metKEGGID;
    targetListKEGG = mapIDsViaMNXref('mets',query1,'KEGG','MetaNetX');
    
    %mapping model.metChEBIID to metanet ID, due to the chebi IDs were saved as
    %CHEBI:##### this format, so before we do the mapping, we need to do some
    %preparation.
    query2 = model.metChEBIID;
    for i = 1:length(query2)
        a = query2{i};
        queryList{i} = a(7:end);
    end
    queryList = transpose(queryList);
    targetListChEBI = mapIDsViaMNXref('mets',queryList,'ChEBI','MetaNetX');
    
    %match to model.mets
    if level == 1
        for i =1:length(model.mets)
            if isempty(model.metMetaNetXID{i})
                if strcmp(targetListKEGG{i,1},targetListChEBI{i,1}) && ~isempty(targetListChEBI{i,1})
                    model.metMetaNetXID{i} = targetListKEGG{i,1};
                end
            end
        end
        
    elseif level == 2
        %load mapping list for yeast7.6MNXmodel from mnx website
        fid  = fopen('../../ComplementaryData/databases/Yeast7.6MNXMetMappingList.tsv','r');
        MNX = textscan(fid,'%s %s','Delimiter','\t','HeaderLines',1);
        MNXmodel.mets  = MNX{1};
        MNXmodel.MetMNXid = MNX{2};
        fclose(fid);
        
        for i =1:length(model.mets)
            for j = 1:length(MNXmodel.mets)
                if isequal(model.mets{i},MNXmodel.mets{j})&& isempty(model.metMetaNetXID{i})
                    model.metMetaNetXID{i} = MNXmodel.MetMNXid{j,1};
                elseif strcmp(targetListKEGG{i,1},targetListChEBI{i,1}) && ~isempty(targetListKEGG{i,1}) && isempty(model.metMetaNetXID{i})
                    model.metMetaNetXID{i} = targetListKEGG{i,1};
                end
            end
        end
    end
    
elseif strcmpi(type,'rxns')
    query3 = model.rxnKEGGID;
    targetListKEGG = mapIDsViaMNXref('rxns',query3,'KEGG','MetaNetX');
    %match to model.mets
    if level == 1
        for i =1:length(model.rxns)
            if isempty(model.rxnMetaNetXID{i})
                model.metMetaNetXID{i} = targetListKEGG{i,1};
            end
        end
        
    elseif level == 2
        %load mapping list for yeast7.6MNXmodel from mnx website
        fid  = fopen('../../ComplementaryData/databases/Yeast7.6MNXRxnMappingList.tsv','r');
        MNX = textscan(fid,'%s %s','Delimiter','\t','HeaderLines',1);
        MNXmodel.rxns  = MNX{1};
        MNXmodel.rxnMNXid = MNX{2};
        fclose(fid);
        
        for i =1:length(model.rxns)
            for j = 1:length(MNXmodel.rxns)
                if isequal(model.rxns{i},MNXmodel.rxns{j})&& isempty(model.rxnMetaNetXID{i})
                    model.rxnMetaNetXID{i} = MNXmodel.rxnMNXid{j,1};
                elseif ~isempty(targetListKEGG{i,1}) && isempty(model.rxnMetaNetXID{i})
                    model.rxnMetaNetXID{i} = targetListKEGG{i,1};
                end
            end
        end
    end
    
else
    EM='Incorrect value of the "type" parameter. Allowed values are "rxns" or "mets"';
    dispEM(EM);
end

%save model
cd ..
saveYeastModel(model)
end

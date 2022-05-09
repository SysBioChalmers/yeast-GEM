function newModel=curateMetsRxnsGenes(model,metsInfo,genesInfo,rxnsCoeffs,rxnsInfo)
% curateMetsRxnsGenes
%   Curate existing and/or add new metabolites, reactions and genes. The
%   information on what metabolites, reactions and/or genes to add are
%   provided in four .tsv files. Templates of these files are given in the
%   data/modelCuration/template folder. Copy these files to your favourite
%   location and modify them to match what you want to curate.
%
%   If the *.tsv files contain metabolites, reactions and/or genes that are
%   already present in the model, then information in the model will be
%   overwritten. Note that this includes empty annotations in the *.tsv
%   files! Metabolites are matched by metaboliteName[comp]; reactions by
%   the stoichiometry of its reactants and products; genes by their gene
%   name. This function can therefore be used to add new entities in the
%   model, or curate those already existing in the model.
%
%   Input:
%       model       RAVEN model structure to be curated
%       metsInfo    relative path to the *.tsv file with metabolite
%                   information. Path to the template file would be:
%                   '../data/modelCuration/template/metsInfo.tsv'. If no
%                   metabolites should be curated, metsInfo should be 'none'.
%       genesInfo   relative path to the *.tsv file with gene
%                   information. Path to the template file would be:
%                   '../data/modelCuration/template/genesInfo.tsv'. If no
%                   genes should be curated, genesInfo should be 'none'.
%       rxnsCoeffs  relative path to the *.tsv file with stoichiometric
%                   coefficients. Path to the template file would be:
%                   '../data/modelCuration/template/rxnsCoeffs.tsv'. If no
%                   rxns should be curated, rxnsCoeffs should be 'none'.
%       rxnsInfo    relative path to the *.tsv file with reaction
%                   information. Path to the template file would be:
%                   '../data/modelCuration/template/rxnsInfo.tsv'. If no
%                   rxns should be curated, rxnsInfo should be 'none'.
%
%   Output:
%       newModel    curated RAVEN model structure

if nargin==4
    error('Provide both a ''rxnsInfo'' and a ''rxnsCoeffs'' file')
end
if nargin<4
    rxnsInfo='none';
    rxnsCoeffs='none';
end
if nargin<3
    genesInfo='none';
end
newModel=model;

%% Metabolites
if ~strcmp(metsInfo,'none')
    fid = fopen(metsInfo);
    raw = textscan(fid,['%q' repmat(' %q',1,17)],'Delimiter','\t');
    fclose(fid);
    metsToAdd.metNames      = raw{1}(2:end);
    metsToAdd.compartments  = raw{2}(2:end);
    metsToAdd.metFormulas   = raw{3}(2:end);
    metsToAdd.metCharges    = cellfun(@str2num,raw{4}(2:end),'UniformOutput',false);
    emptyEntry = cellfun(@isempty,metsToAdd.metCharges);
    if all(emptyEntry)
        metsToAdd           = rmfield(metsToAdd,'metCharges');
    elseif any(emptyEntry) % If some charges are given, assume 0 for those without charges specified
        metsToAdd.metCharges(emptyEntry) = {0};
        metsToAdd.metCharges = cell2mat(metsToAdd.metCharges);
    else
        metsToAdd.metCharges = cell2mat(metsToAdd.metCharges);
    end
    metsToAdd.inchis        = raw{5}(2:end);
    metsToAdd.metNotes      = raw{6}(2:end);
    
    % Check if metabolite already exists (check by metName[comp])
    existingMets            = [];
    existingMetsIdx         = [];
    newMets                 = [];
    for i=1:numel(metsToAdd.metNames)
        metIdx              = find(strcmp(newModel.metNames,metsToAdd.metNames{i}));
        existMet            = strcmp(metsToAdd.compartments{i},newModel.comps(newModel.metComps(metIdx)));
        if any(existMet)
            existingMets	= [existingMets, i];
            existingMetsIdx = [existingMetsIdx, metIdx(existMet)];
        else
            newMets	= [newMets, i];
        end
    end
    
    % Overwrite annotation of existing entries
    if any(existingMets)
        if isfield(newModel,'metFormulas')
            newModel.metFormulas(existingMetsIdx)  = metsToAdd.metFormulas(existingMets);
        end
        if isfield(newModel,'metCharges') && isfield(metsToAdd,'metCharges')
            newModel.metCharges(existingMetsIdx)   = metsToAdd.metCharges(existingMets);
        end
        if isfield(newModel,'inchis')
            newModel.inchis(existingMetsIdx)       = metsToAdd.inchis(existingMets);
        end
        newModel = extractAndAddMiriam(newModel,raw(7:end),existingMets,existingMetsIdx,'met');
        warning(['The following metabolites are already present in the model, '...
            'their annotation will be overwritten to match the metsInfo file. '...
            'If you do not particularly want to curate their annotations, it '...
            'would be better to removes these metabolites from metsInfo:\n\t\t%s'],...
            strjoin(strcat(metsToAdd.metNames(existingMets),'[',metsToAdd.compartments(existingMets),']'),'\n\t\t'));
    end
    
    % Continue with new metabolites
    if numel(newMets)>0
        metsToAdd.metNames(existingMets)        = [];
        metsToAdd.compartments(existingMets)    = [];
        metsToAdd.metFormulas(existingMets)     = [];
        if all(cellfun(@isempty,metsToAdd.metFormulas))
            metsToAdd = rmfield(metsToAdd,'metFormulas');
        end
        if isfield(metsToAdd,'metCharges')
            metsToAdd.metCharges(existingMets)	= [];
        end
        metsToAdd.inchis(existingMets)          = [];
        if all(cellfun(@isempty,metsToAdd.inchis))
            metsToAdd = rmfield(metsToAdd,'inchis');
        end
        metsToAdd.metNotes(existingMets)        = [];
        if all(cellfun(@isempty,metsToAdd.metNotes))
            metsToAdd = rmfield(metsToAdd,'metNotes');
        end
        % Add metabolites
        newModel        = addMets(newModel,metsToAdd,true,'s_');
        addedIdx        = numel(newModel.mets)-numel(newMets)+1:numel(newModel.mets);
        newModel        = extractAndAddMiriam(newModel,raw(7:end),newMets,addedIdx,'met');
    end
end
%% Genes
if ~strcmp(genesInfo,'none')
    % Gather all data, first about reaction stoichiometries...
    fid = fopen(genesInfo);
    raw = textscan(fid,['%q' repmat(' %q',1,10)],'Delimiter','\t');
    fclose(fid);
    genesToAdd.genes            = raw{1}(2:end);
    genesToAdd.geneShortNames   = raw{2}(2:end);

    [~,notNewGene,existingGene] = intersect(genesToAdd.genes,newModel.genes);
    if ~isempty(notNewGene)
        warning(['The following genes are already present in the model, their annotation '...
            'will be overwritten to match the genesInfo file: \n\t\t%s'],...
            strjoin(genesToAdd.genes(notNewGene),'\n\t\t'));
        if isfield(newModel,'geneShortNames')
            newModel.geneShortNames(existingGene) = genesToAdd.geneShortNames(notNewGene);
        end
        newModel = extractAndAddMiriam(newModel,raw(3:end),notNewGene,existingGene,'gene');
    end
    
    % Continue with new genes
    toAdd               = 1:numel(genesToAdd.genes);
    toAdd(notNewGene)   = [];
    if ~isempty(toAdd)
        genesToAdd.genes            = genesToAdd.genes(toAdd);
        genesToAdd.geneShortNames   = genesToAdd.geneShortNames(toAdd);
        if all(cellfun(@isempty,genesToAdd.geneShortNames))
            genesToAdd = rmfield(genesToAdd,'geneShortNames');
        end
        % Add genes
        newModel = addGenesRaven(newModel,genesToAdd);
        addedIdx = numel(newModel.genes)-numel(toAdd)+1:numel(newModel.genes);
        newModel = extractAndAddMiriam(newModel,raw(3:end),toAdd,addedIdx,'gene');
    end
end
%% Reactions
if ~any(strcmp({rxnsCoeffs,rxnsInfo},'none'))
    % Gather all data, first about reaction stoichiometries...
    fid = fopen(rxnsCoeffs);
    raw = textscan(fid,'%q %q %q %q %f','Delimiter','\t','HeaderLines',1);
    fclose(fid);
    rxnCheck.coeffsIdx = str2double(raw{1});
    rxns.rxnNames   = raw{2};
    rxns.metNames   = raw{3};
    rxns.comps      = raw{4};
    rxns.coeffs     = raw{5};
    rxnCheck.coeffs = strcat(raw{1},'***',rxns.rxnNames);

    % ... and then additional reaction-specific data.
    fid = fopen(rxnsInfo);
    raw = textscan(fid,['%q' repmat(' %q',1,20)],'Delimiter','\t');
    fclose(fid);
    rxnCheck.rxnsIdx = str2double(raw{1}(2:end));
    rxnCheck.rxns = strcat(raw{1}(2:end),'***',raw{2}(2:end));
    
    notMatching = setxor(rxnCheck.coeffs,rxnCheck.rxns);
    if numel(notMatching)>1
        error(['The following reactions ánd/or their indices are not matched '...
            'between the rxnsInfo and rxnsCoeffs files:\n\t\t%s'],...
            strjoin(regexprep(notMatching,'^\d+***',''),'\n\t\t'))
    end
    
    rxnsToAdd.rxnNames              = raw{2}(2:end);
    rxnsToAdd.grRules               = raw{3}(2:end);
    rxnsToAdd.lb                    = cellfun(@str2num,raw{4}(2:end));
    rxnsToAdd.ub                    = cellfun(@str2num,raw{5}(2:end));
    rxnsToAdd.rev                   = cellfun(@str2num,raw{6}(2:end));
    rxnsToAdd.subSystems            = raw{7}(2:end);
    rxnsToAdd.eccodes               = raw{8}(2:end);
    rxnsToAdd.rxnNotes              = raw{9}(2:end);
    rxnsToAdd.rxnReferences         = raw{10}(2:end);
    rxnsToAdd.rxnConfidenceScores   = cellfun(@str2num,raw{11}(2:end),'UniformOutput',false);
    emptyEntry = cellfun(@isempty,rxnsToAdd.rxnConfidenceScores);
    if all(emptyEntry)
        rxnsToAdd           = rmfield(rxnsToAdd,'rxnConfidenceScores');
    elseif any(emptyEntry) % If some rxnConfidenceScores are given, assume 0 for those without rxnConfidenceScores specified
        rxnsToAdd.rxnConfidenceScores(emptyEntry) = {0};
        rxnsToAdd.rxnConfidenceScores = cell2mat(rxnsToAdd.rxnConfidenceScores);
    else
        rxnsToAdd.rxnConfidenceScores = cell2mat(rxnsToAdd.rxnConfidenceScores);
    end
    
    existingRxn=[];
    notNewRxn=[];
    for i=1:numel(rxnCheck.rxnsIdx)
        rxnRows                     = find(rxnCheck.rxnsIdx(i)==rxnCheck.coeffsIdx);
        rxnsToAdd.mets{i,1}         = cell(1,1);
        rxnsToAdd.stoichCoeffs{i,1} = cell(1,1);
        for j=1:numel(rxnRows)
            newMetComp = [rxns.metNames{rxnRows(j)}, '[', rxns.comps{rxnRows(j)}, ']'];
            try
                newMetComp = getIndexes(newModel,newMetComp,'metcomps');
            catch
                error(['Not all metabolites in reaction "%s" are present in the '...
                    'model or in the provided table with new metabolites.'],rxnsToAdd.rxnNames{i})
            end
            rxnsToAdd.mets{i}(j) = newModel.mets(newMetComp);
            rxnsToAdd.stoichCoeffs{i}{j} = rxns.coeffs(rxnRows(j));
        end
        rxnsToAdd.stoichCoeffs{i}=cell2mat(rxnsToAdd.stoichCoeffs{i});
        %Check if the reaction not already exists, by stoichiometry of its
        %products and reactants
        modelCoeffs         = transpose(newModel.S(getIndexes(newModel,rxnsToAdd.mets{i},'mets'),:));
        modelCoeffs2        = find(~all(modelCoeffs==0,2));
        modelCoeffs         = modelCoeffs(modelCoeffs2,:);
        [~, duplicateRxn]   = intersect(modelCoeffs,rxnsToAdd.stoichCoeffs{i},'rows');
        if ~isempty(duplicateRxn)
            existingRxn     = [existingRxn,modelCoeffs2(duplicateRxn)];
            notNewRxn       = [notNewRxn,i];
        end
    end
    
    % If some reactions already existed, then replace its information
    if ~isempty(existingRxn)
        warning(['The following reactions are the same as existing model '...
            'reactions, their annotation will be overwritten to match the '...
            'rxnCoeffs file: \n\t\t(Existing reaction ID): New reaction name\n\t\t%s'],...
            strjoin(strcat('(',newModel.rxns(existingRxn),{'): '},rxnsToAdd.rxnNames(notNewRxn)),'\n\t\t'));
        newModel.rxnNames(existingRxn) = rxnsToAdd.rxnNames(notNewRxn);
        newModel.lb(existingRxn)       = rxnsToAdd.lb(notNewRxn);
        newModel.ub(existingRxn)       = rxnsToAdd.ub(notNewRxn);
        newModel.rev(existingRxn)      = rxnsToAdd.rev(notNewRxn);
        if isfield(newModel,'subSystems')
            newModel.subSystems(existingRxn)    = rxnsToAdd.subSystems(notNewRxn);
        end
        if isfield(newModel,'eccodes')
            newModel.eccodes(existingRxn)       = rxnsToAdd.eccodes(notNewRxn);
        end
        if isfield(newModel,'rxnNotes')
            newModel.rxnNotes(existingRxn)      = rxnsToAdd.rxnNotes(notNewRxn);
        end
        if isfield(newModel,'rxnConfidenceScores') && isfield(rxnsToAdd,'rxnConfidenceScores')
            newModel.rxnConfidenceScores(existingRxn) = rxnsToAdd.rxnConfidenceScores(notNewRxn);
        end
        emptyEntries = cellfun(@isempty,rxnsToAdd.grRules);
        if ~all(emptyEntries)
            newModel    = changeGrRules(newModel,newModel.rxns(existingRxn(~emptyEntries)),rxnsToAdd.grRules(notNewRxn(~emptyEntries)),true);
        end        
        newModel    = extractAndAddMiriam(newModel,raw(11:end),notNewRxn,existingRxn,'rxn');
    end
    
    % Continue with new reactions
    toAdd               = 1:numel(rxnsToAdd.rxnNames);
    toAdd(notNewRxn)    = [];
    if ~isempty(toAdd)
        rxnsToAdd.rxnNames      = rxnsToAdd.rxnNames(toAdd);
        rxnsToAdd.lb            = rxnsToAdd.lb(toAdd);
        rxnsToAdd.ub            = rxnsToAdd.ub(toAdd);
        rxnsToAdd.rev           = rxnsToAdd.rev(toAdd);
        rxnsToAdd.subSystems    = rxnsToAdd.subSystems(toAdd);
        if all(cellfun(@isempty,rxnsToAdd.subSystems))
            rxnsToAdd = rmfield(rxnsToAdd,'subSystems');
        end
        rxnsToAdd.eccodes       = rxnsToAdd.eccodes(toAdd);
        if all(cellfun(@isempty,rxnsToAdd.eccodes))
            rxnsToAdd = rmfield(rxnsToAdd,'eccodes');
        end
        rxnsToAdd.rxnNotes      = rxnsToAdd.rxnNotes(toAdd);
        if all(cellfun(@isempty,rxnsToAdd.rxnNotes))
            rxnsToAdd = rmfield(rxnsToAdd,'rxnNotes');
        end
        if isfield(rxnsToAdd,'rxnConfidenceScores')
            rxnsToAdd.rxnConfidenceScores = rxnsToAdd.rxnConfidenceScores(toAdd);
        end
        rxnsToAdd.rxnReferences = rxnsToAdd.rxnReferences(toAdd);
        if all(cellfun(@isempty,rxnsToAdd.rxnReferences))
            rxnsToAdd = rmfield(rxnsToAdd,'rxnReferences');
        end
        rxnsToAdd.grRules       = rxnsToAdd.grRules(toAdd);
        rxnsToAdd.mets          = rxnsToAdd.mets(toAdd);
        rxnsToAdd.stoichCoeffs  = rxnsToAdd.stoichCoeffs(toAdd);
        rxnsToAdd.rxns          = generateNewIds(model,'rxns','r_',numel(rxnsToAdd.rxnNames));
        
        newModel        = addRxns(newModel,rxnsToAdd,1,[],false,false);
        rxnsModelIdx    = numel(newModel.rxns)-numel(toAdd)+1:numel(newModel.rxns);
        newModel        = extractAndAddMiriam(newModel,raw(11:end),toAdd,rxnsModelIdx,'rxn');
    end
end
end

function newModel = extractAndAddMiriam(model,raw,inputIndex,modelIndex,type)
newModel=model;
miriamName              = cell(numel(inputIndex),1);
miriamValues            = cell(1,numel(inputIndex));
for i=1:numel(raw)
    miriamName{i}                       = raw{i}{1};
    miriamValues(1:numel(inputIndex),i) = raw{i}(inputIndex+1);
end
emptyMiriam     =   all(cellfun(@isempty,miriamValues),1);
miriamName(emptyMiriam)     = [];
miriamValues(:,emptyMiriam) = [];
if ~isfield(newModel,[type 'Miriams']);
    newModel.([type 'Miriams'])=cell(numel(newModel.([type 's'])),1);
end
if ~isempty(miriamName)
    for i=1:numel(miriamName)
        newModel = editMiriam(newModel,type,modelIndex,miriamName{i},miriamValues(:,i),'replace');
    end
end
end

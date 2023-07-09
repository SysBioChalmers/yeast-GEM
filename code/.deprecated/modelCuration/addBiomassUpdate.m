function model = addBiomassUpdate(model,data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addBiomassUpdate
% Adds rxns and metabolites for BiomassUpdate.
% Input: model, Biomass_newRxnMatrix.tsv, Biomass_newRxnProp.tsv,
%        Biomass_newRxnMet.tsv.
%
% NOTE: changeGeneAssociation.m is a function from COBRA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Load stoichiometry data:
cd ..
fid = fopen('../data/modelCuration/Biomass_newRxnMatrix.tsv');
newreaction = textscan(fid,'%s %s %s %s %s','Delimiter','\t','HeaderLines',1);
matrix.rxnIDs  = newreaction{1};
matrix.metcoef = cellfun(@str2num, newreaction{2});
matrix.metIDs  = newreaction{3};
matrix.mettype = newreaction{4};
matrix.metcompartments = newreaction{5};
fclose(fid);

% Load rxn properties data:
fid  = fopen('../data/modelCuration/Biomass_newRxnProp.tsv','r');
rev = textscan(fid,'%s %s %s %s %s %s %s','Delimiter','\t','HeaderLines',1);
newrxn.ID  = rev{1};
newrxn.Rev = cellfun(@str2num, rev{2});
newrxn.GPR = rev{3};
newrxn.rxnNames      = rev{4};
newrxn.rxnECNumbers  = rev{5};
newrxn.rxnKEGGID     = rev{6};
newrxn.rxnNotes      = rev{7};
newrxn.rxnMetaNetXID = newrxn.ID;
for i = 1:length(newrxn.rxnMetaNetXID)
    if ~startsWith(newrxn.rxnMetaNetXID{i},'MNXR')
        newrxn.rxnMetaNetXID{i} = '';
    end
end
newrxn.rxnMetaNetXID = regexprep(newrxn.rxnMetaNetXID,'_cv','');
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
[m, ~]      = size(CONValldata);
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
fid = fopen('../../data/modelCuration/Biomass_newRxnMet.tsv');
newmet_annot         = textscan(fid,'%s %s %s %s %s %s %s','Delimiter','\t','HeaderLines',1);
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
        model.metNotes{metID}      = 'NOTES:added for BiomassUpdate (PR #174)';
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
        'geneRule',newrxn.GPR{i},...
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
    model.rxnNotes{rxnIndex} = 'NOTES:added for BiomassUpdate (PR #174)';
    if ~isempty(model.rules{rxnIndex})
        model.rxnConfidenceScores(rxnIndex) = 2;   %reactions has GPR rules
    elseif contains(model.rxnNames(rxnIndex),'exchange')
        model.rxnConfidenceScores(rxnIndex) = NaN; %exchange rxns
    else
        model.rxnConfidenceScores(rxnIndex) = 1;   %reactions without gene but needed for modelling
    end
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Introduce two PseudoRxn for cofactors and ions
group = unique(data.groups);
for i = 1:length(group)
%get new metID for pseudo mets :ion and cofactor
    cd ../otherChanges
    metname = [group{i} ' [cytoplasm]'];
    [~,index] = ismember(metname,model.metNames);
    if index ~= 0
        metID = model.mets(index);
    elseif index == 0
        metID  = ['s_' getNewIndex(model.mets) '[c]'];
        model = addMetabolite(model,metID, ...
            'metName',metname);
    end
    [~,index] = ismember(metID,model.mets);
    model.metNotes{index}      = 'NOTES:added for BiomassUpdate (PR #174)';

% Get coefficients for pseudo rxns: ion pseudoreaction, cofactor pseudoreaction
    newps.rxnName = {[group{i},' pseudoreaction']};
    index = strcmp(group{i},data.groups);
    newps.metIDs = data.mets(index);
    newps.metcoef = data.abundances(index)*-1;

% introduce two pseudo mets ion [cytoplasm], cofactor [cytoplasm] into the reaction
    newps.metIDs = [newps.metIDs;metID];
    newps.metcoef = [newps.metcoef;1];
    [~,index] = ismember(newps.rxnName,model.rxnNames);
    if index == 0
        rxnID = ['r_' getNewIndex(model.rxns)];
    else
        rxnID = model.rxns{index};
    end
    [model,rxnIndex] = addReaction(model, rxnID,...
                       'reactionName', char(newps.rxnName),...
                       'metaboliteList',newps.metIDs,...
                       'stoichCoeffList',newps.metcoef,...
                       'reversible',false,...
                       'checkDuplicate',1);
     cd ../modelCuration/
    [EnergyResults,RedoxResults] = CheckEnergyProduction(model,{['r_' newID]},EnergyResults,RedoxResults);
    [MassChargeresults] = CheckBalanceforSce(model,{['r_' newID]},MassChargeresults);
    if isempty(rxnIndex)
        rxnIndex = strcmp(model.rxns,rxnID);
    end
    model.rxnNotes{rxnIndex} = 'NOTES:added for BiomassUpdate (PR #174)';
    model.rxnConfidenceScores(rxnIndex) = 1; %needed for modeling
    model = rmfield(model,'grRules');
end

% Change original biomass equation
UnusedMets = {'sulphate [cytoplasm]','riboflavin [cytoplasm]', 'heme a [cytoplasm]'};
[~,UnusedMets_index] = ismember(UnusedMets,model.metNames);
AddMets = {'cofactor [cytoplasm]','ion [cytoplasm]'};
[~,AddMets_index] = ismember(AddMets,model.metNames);
biomassRxn = 'r_4041';
[~,biomassRxn_index] = ismember(biomassRxn,model.rxns);
model.S(UnusedMets_index,biomassRxn_index) = 0;
model.S(AddMets_index,biomassRxn_index) = -1;

% Save model:
if isfield(model,'grRules')
model = rmfield(model,'grRules');
end
model = minimal_Y6(model);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script is used to sort reactions by using deltaG from various 
% sources, based on the defined deltaG threshold of ±30 kJ/mol.
% Furthermore, reversibility index is calculated using get_lnRevIdx.m and
% reaction direction information is compiled
%
% Note: script is added for documentation purposes, since the output
%       rxnDirectionInfo.tsv has already been edited and added
%       into into repository
%
% Input: model, grpContribution_deltaG.tsv, MetaCyc_deltaG.tsv,
%        rev_index.mat
% Note: grpContribution_deltaG.tsv is taken from supplmentary data in
%       https://doi.org/10.1529/biophysj.107.124784
%
% Output: allrxns (also generates rxnDirectionInfo.tsv if the last segment 
%                  of the script is executed)     
%

%Load model
cd ..
model = loadYeastModel;

%List of reactions which involves nucleotides
NTP = {'ATP';'GTP';'CTP';'UTP';'TTP';'ADP';'GDP';'CDP';'UDP';'TDP';'AMP';'GMP';'CMP';'UMP';'TMP';...
    'dATP';'dGTP';'dCTP';'dUTP';'dTTP';'dADP';'dGDP';'dCDP';'dUDP';'dTDP';'dAMP';'dGMP';'dCMP';'dUMP';'dTMP'};

modelR = ravenCobraWrapper(model);
met_idx = find(ismember(modelR.metNames,NTP));
tmp = model.S(met_idx,:);
exclrxnList = model.rxns(any(tmp,1));
printRxnFormula(model,'rxnAbbrList',exclrxnList,'metNameFlag',true)
excl_rxns = exclrxnList;

%List of exchange reactions
exchangeRxns = findExcRxns(model);
excl_rxns = [excl_rxns;model.rxns(exchangeRxns)];
excl_rxns = unique(excl_rxns);

%Generate list of reactions to check
checkRxns = setdiff(model.rxns,excl_rxns);

%Categorize reactions via deltaG from group contribution method
fid2 = fopen('../data/modelCuration/grpContribution_deltaG.tsv');
format = repmat('%s ',1,4);
format = strtrim(format);
rxn_temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',0);
for i = 1:length(rxn_temp)
    grpContribution(:,i) = rxn_temp{i};
end
commentLines = startsWith(grpContribution(:,1),'#');
grpContribution(commentLines,:) = []; % rxnKEGGID deltaG Uncertainty Notes
fclose(fid2);

KEGGID = grpContribution(:,1);
deltaG_grpCont = grpContribution(:,2);
uncertainty_grpCont = grpContribution(:,3);
notes = grpContribution(:,4);

bwd_grpCont = cell(length(model.rxns),4);
fwd_grpCont = cell(length(model.rxns),4);
rev_grpCont = cell(length(model.rxns),4);
no_deltaG = cell(length(model.rxns),4);
noMapping = cell(length(model.rxns),3);

[~,rxn_idx] = ismember(checkRxns,model.rxns);
rxnKEGGID = model.rxnKEGGID(rxn_idx);

cd modelCuration/
for i = 1:length(rxnKEGGID)
    if ~ismissing(rxnKEGGID(i)) && contains(rxnKEGGID(i),'R')
        idx = find(ismember(KEGGID,rxnKEGGID(i)));
        if isempty(idx)
            %rxnKEGGID not found in grpContribution_deltaG.tsv
            arrayidx = find(cellfun('isempty',noMapping),1);
            noMapping(arrayidx,:) = [checkRxns(i),rxnKEGGID(i),'rxnKEGGID not found in grpContribution_deltaG.tsv'];
        elseif ~ismember(deltaG_grpCont(idx),'Not calculated') && ~ismember(uncertainty_grpCont(idx),'Not calculated')
            %deltaG range for rev_GrpKEGG reaction: -30 <= mod_deltaG <= 30
            deltaG = str2num(char(deltaG_grpCont(idx)));
            lnRevIdx = get_lnRevIdx(deltaG,checkRxns(i),model);
            uncertainty = str2num(char(uncertainty_grpCont(idx)));
            if deltaG > 0
                mod_deltaG = deltaG - uncertainty;
                if mod_deltaG > 30
                    arrayidx = find(cellfun('isempty',bwd_grpCont),1);
                    bwd_grpCont(arrayidx,:) = [checkRxns(i),deltaG,uncertainty,lnRevIdx];
                else
                    arrayidx = find(cellfun('isempty',rev_grpCont),1);
                    rev_grpCont(arrayidx,:) = [checkRxns(i),deltaG,uncertainty,lnRevIdx];
                end
            elseif deltaG < 0
                mod_deltaG = deltaG + uncertainty;
                if mod_deltaG < -30
                    arrayidx = find(cellfun('isempty',fwd_grpCont),1);
                    fwd_grpCont(arrayidx,:) = [checkRxns(i),deltaG,uncertainty,lnRevIdx];
                else
                    arrayidx = find(cellfun('isempty',rev_grpCont),1);
                    rev_grpCont(arrayidx,:) = [checkRxns(i),deltaG,uncertainty,lnRevIdx];
                end
            else
                arrayidx = find(cellfun('isempty',rev_grpCont),1);
                rev_grpCont(arrayidx,:) = [checkRxns(i),deltaG,uncertainty,lnRevIdx];
            end
        else %deltaG/uncertainty not calculated
            arrayidx = find(cellfun('isempty',no_deltaG),1);
            no_deltaG(arrayidx,:) = [checkRxns(i),'Not calculated','-',notes(idx)];
        end
    else %no rxnKEGGID in model
        arrayidx = find(cellfun('isempty',noMapping),1);
        noMapping(arrayidx,:) = [checkRxns(i),'No rxnKEGGID in model','-'];
    end
end

%remove empty cells
empties = find(cellfun('isempty',bwd_grpCont(:,1)));
bwd_grpCont(empties,:) = [];
empties = find(cellfun('isempty',fwd_grpCont(:,1)));
fwd_grpCont(empties,:) = [];
empties = find(cellfun('isempty',rev_grpCont(:,1)));
rev_grpCont(empties,:) = [];
empties = find(cellfun('isempty',noMapping(:,1)));
noMapping(empties,:) = [];
empties = find(cellfun('isempty',no_deltaG(:,1)));
no_deltaG(empties,:) = [];

%Categorize reactions via deltaG from MetaCyc database
cd ../missingFields/
[~,rxn_idx] = ismember(checkRxns,model.rxns);
rxnMetaNetXID = model.rxnMetaNetXID(rxn_idx);
xref_metacyc = mapIDsViaMNXref('rxns',rxnMetaNetXID,'MetaNetX','MetaCyc');
xref_metacyc(:,2) = checkRxns; %match rxnMetaCycID with checkRxns
empties = find(cellfun('isempty',xref_metacyc(:,1)));
idx_check2 = empties;
xref_metacyc(empties,:) = [];

%load metaCycRxns.mat via getRxnsFromMetaCyc from RAVEN directory
MetaCyc_rxnInfo = [];
try
    MetaCyc_rxnInfo = getRxnsFromMetaCyc;
catch
    warning('RAVEN repository is not cloned or added to MATLAB path, unable to use function getRxnsFromMetaCyc');
end

%load MetaCyc_deltaG.tsv as getRxnsFromMetaCyc does not include deltaG data
cd ..
fid2 = fopen('../data/modelCuration/MetaCyc_deltaG.tsv');
format = repmat('%s ',1,2);
format = strtrim(format);
rxn_temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',0);
for i = 1:length(rxn_temp)
    MetaCyc_deltaG(:,i) = rxn_temp{i};
end
commentLines = startsWith(MetaCyc_deltaG(:,1),'#');
MetaCyc_deltaG(commentLines,:) = []; % rxnMetaCycID deltaG
fclose(fid2);
[~,rxnIdx2] = ismember(xref_metacyc(:,1), MetaCyc_deltaG(:,1));
[~,rxnIdx3] = ismember(xref_metacyc(:,1), MetaCyc_rxnInfo.rxns);

%preallocate cell arrays
bwd_MetaCyc{length(model.rxns),4} = [];
fwd_MetaCyc{length(model.rxns),4} = [];
rev_MetaCyc{length(model.rxns),4} = [];
no_deltaG2{length(model.rxns),4} = [];
noMapping2{length(model.rxns),4} = [];

cd modelCuration/
if ~isempty(MetaCyc_rxnInfo) %check if MetaCyc_rxnInfo is generated successfully
    for i = 1:length(rxnIdx2)
        if rxnIdx2(i) ~= 0
            if ~isempty(MetaCyc_deltaG{rxnIdx2(i),2})
                deltaG = str2num(char(MetaCyc_deltaG(rxnIdx2(i),2))); %MetaCyc reaction.dat only provided Gibbs-0 information (no uncertainty given)
                lnRevIdx = get_lnRevIdx(deltaG,xref_metacyc(i,2),model);
                if deltaG > 30
                    arrayidx = find(cellfun('isempty',bwd_MetaCyc),1);
                    if rxnIdx3(i) ~= 0
                        if MetaCyc_rxnInfo.rev(rxnIdx3(i)) == 1
                            rxnDirection = {'='};
                        elseif MetaCyc_rxnInfo.rev(rxnIdx3(i)) == 0
                            rxnDirection = {'>'};
                        end
                        bwd_MetaCyc(arrayidx,:) = [xref_metacyc(i,2),deltaG,lnRevIdx,rxnDirection];
                    else
                        bwd_MetaCyc(arrayidx,:) = [xref_metacyc(i,2),deltaG,lnRevIdx,'-'];
                    end
                elseif deltaG < -30
                    arrayidx = find(cellfun('isempty',fwd_MetaCyc),1);
                    if rxnIdx3(i) ~= 0
                        if MetaCyc_rxnInfo.rev(rxnIdx3(i)) == 1
                            rxnDirection = {'='};
                        elseif MetaCyc_rxnInfo.rev(rxnIdx3(i)) == 0
                            rxnDirection = {'>'};
                        end
                        fwd_MetaCyc(arrayidx,:) = [xref_metacyc(i,2),deltaG,lnRevIdx,rxnDirection];
                    else
                        fwd_MetaCyc(arrayidx,:) = [xref_metacyc(i,2),deltaG,lnRevIdx,'-'];
                    end
                else
                    arrayidx = find(cellfun('isempty',rev_MetaCyc),1);
                    if rxnIdx3(i) ~= 0
                        if MetaCyc_rxnInfo.rev(rxnIdx3(i)) == 1
                            rxnDirection = {'='};
                        elseif MetaCyc_rxnInfo.rev(rxnIdx3(i)) == 0
                            rxnDirection = {'>'};
                        end
                        if isempty(deltaG)
                            deltaG = {[]};
                            lnRevIdx = {[]};
                        end
                        rev_MetaCyc(arrayidx,:) = [xref_metacyc(i,2),deltaG,lnRevIdx,rxnDirection];
                    else
                        if isempty(deltaG)
                            deltaG = {[]};
                            lnRevIdx = {[]};
                        end
                        rev_MetaCyc(arrayidx,:) = [xref_metacyc(i,2),deltaG,lnRevIdx,{''}];
                    end
                end
            else
                if rxnIdx3(i) ~= 0 %if no deltaG, sort based on rxnDirection (if available)
                    if MetaCyc_rxnInfo.rev(rxnIdx3(i)) == 1
                        rxnDirection = {'='};
                        arrayidx = find(cellfun('isempty',rev_MetaCyc),1);
                        rev_MetaCyc(arrayidx,:) = [xref_metacyc(i,2),{[]},{[]},rxnDirection];
                    elseif MetaCyc_rxnInfo.rev(rxnIdx3(i)) == 0
                        rxnDirection = {'>'};
                        arrayidx = find(cellfun('isempty',fwd_MetaCyc),1);
                        fwd_MetaCyc(arrayidx,:) = [xref_metacyc(i,2),{[]},{[]},rxnDirection];
                    end
                else
                    arrayidx = find(cellfun('isempty',no_deltaG2),1);
                    no_deltaG2(arrayidx,:) = [xref_metacyc(i,2),xref_metacyc(i,1),'No deltaG in MetaCyc database','-'];
                end
            end
        else
            arrayidx = find(cellfun('isempty',noMapping2),1);
            noMapping2(arrayidx,:) = [xref_metacyc(i,2),xref_metacyc(i,1),'MetaCycID not found in MetaCyc database','-'];
        end
    end
    
    for i = 1:length(idx_check2)
        arrayidx = find(cellfun('isempty',noMapping2),1);
        noMapping2(arrayidx,:) = [model.rxns(idx_check2(i)),'-','No MNXRID mapped to rxnMetaCycID','-'];
    end
    
    %remove empty cells
    empties = find(cellfun('isempty',bwd_MetaCyc(:,1)));
    bwd_MetaCyc(empties,:) = [];
    empties = find(cellfun('isempty',fwd_MetaCyc(:,1)));
    fwd_MetaCyc(empties,:) = [];
    empties = find(cellfun('isempty',rev_MetaCyc(:,1)));
    rev_MetaCyc(empties,:) = [];
    empties = find(cellfun('isempty',no_deltaG2(:,1)));
    no_deltaG2(empties,:) = [];
    empties = find(cellfun('isempty',noMapping2(:,1)));
    noMapping2(empties,:) = [];
else
    warning('Checking of reaction direction with MetaCyc IDs unsuccessful as MetaCyc_rxnInfo is not found');
end
    
%Categorize reactions via deltaG from modelSEED database
cd ../missingFields/
[~,rxn_idx] = ismember(checkRxns,model.rxns);
rxnMetaNetXID = model.rxnMetaNetXID(rxn_idx);
xref_seed = mapIDsViaMNXref('rxns', rxnMetaNetXID,'MetaNetX','SEED');
xref_seed(:,2) = checkRxns; %match rxnSEEDID with checkRxns
empties = find(cellfun('isempty',xref_seed(:,1)));
idx_check3 = empties;
xref_seed(empties,:) = [];

%download reactions.tsv from GitHub repository ModelSEED/ModelSEEDDatabase
cd ../modelCuration/
try
    websave('reactions.tsv','https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/dev/Biochemistry/reactions.tsv');
catch
    warning('reactions.tsv was not successfully downloaded, check if directory for reactions.tsv has changed on github.com/ModelSEED/ModelSEEDDatabase/Biochemistry');
end
fileDir = dir('reactions.tsv');

if ~isempty(fileDir) %check if reactions.tsv is present in current directory   
    fid2 = fopen('reactions.tsv');
    format = repmat('%s ',1,22);
    format = strtrim(format);
    rxn_temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',0);
    for i = 1:length(rxn_temp)
        seed_rxnInfo(:,i) = rxn_temp{i};
    end
    seed_rxnInfo(1,:) = []; %remove header in 1st line
    fclose(fid2);
    
    rxnSEEDID = seed_rxnInfo(:,1);
    rxnDefinition = seed_rxnInfo(:,7);
    rxnDirection = seed_rxnInfo(:,10);
    
    %preallocate cell arrays
    bwd_SEED{length(model.rxns),5} = [];
    fwd_SEED{length(model.rxns),5} = [];
    rev_SEED{length(model.rxns),5} = [];
    no_deltaG3{length(model.rxns),5} = [];
    noMapping3{length(model.rxns),3} = [];
    
    [~,rxnIdx4] = ismember(xref_seed(:,1), rxnSEEDID);
    empties3 = find(rxnIdx4==0);
    for i = 1:length(empties3)
        noMapping3(i,:) = [xref_seed(empties3(i),2),xref_seed(empties3(i),1),'No matching of rxnSEEDID in modelSEED database'];
    end
    xref_seed(empties3,:) = []; %only retain rxnSEEDID found in model and contains deltaG in modelSEED databse
    rxnIdx4(empties3,:) = [];
    
    cd ../missingFields/
    xref_seedmets = mapIDsViaMNXref('mets', model.metMetaNetXID,'MetaNetX','SEED');
    cd ../modelCuration/
    for i = 1:length(rxnIdx4)
        if ~isempty(seed_rxnInfo{rxnIdx4(i),15})
            deltaG = str2num(char(seed_rxnInfo(rxnIdx4(i),15)));
            lnRevIdx = get_lnRevIdx(deltaG,xref_seed(i,2),model);
            error = str2num(char(seed_rxnInfo(rxnIdx4(i),16)));
            %check if metabolites are on the same side in both model and SEED
            %then make changes to rxnDirection if necessary
            if strcmp(rxnDirection(rxnIdx4(i)),'<') || strcmp(rxnDirection(rxnIdx4(i)),'>')
                SEED_rxnformula = strrep(rxnDefinition(rxnIdx4(i)),'(','');
                SEED_rxnformula = strrep(SEED_rxnformula,')','');
                SEED_rxnformula = strrep(SEED_rxnformula,' <= ',' <=> '); %make all backward reactions reversible for constructS to work
                [SEED_coef, SEED_mets] = constructS(SEED_rxnformula);
                mets_temp = split(SEED_mets,'[');
                SEED_mets = mets_temp(:,1);
                
                idx_temp = find(ismember(model.rxns,xref_seed(i,2)));
                coef = model.S(:,idx_temp);
                mets_temp = (find(model.S(:,idx_temp)));
                model_mets = xref_seedmets(mets_temp);
                model_coef = coef(mets_temp);
                if length(SEED_mets) > length(model_mets)
                    [idx,idx2] = ismember(model_mets,SEED_mets);
                    if all(idx2)
                        SEED_coef = SEED_coef(idx2);
                        if isequal(SEED_coef,model_coef)
                            %no change to rxnDirection
                        elseif isequal(sign(SEED_coef.*(-1)),sign(model_coef)) && strcmp(rxnDirection(rxnIdx4(i)),'>')
                            rxnDirection(rxnIdx4(i)) = {'<'};
                            def_temp = split(rxnDefinition(rxnIdx4(i)), ' => ');
                            rxnDefinition(rxnIdx4(i)) = strcat(def_temp(3),def_temp(2),def_temp(1));
                        elseif isequal(sign(SEED_coef.*(-1)),sign(model_coef)) && strcmp(rxnDirection(rxnIdx4(i)),'<')
                            rxnDirection(rxnIdx4(i)) = {'>'};
                            def_temp = split(rxnDefinition(rxnIdx4(i)), ' <= ');
                            rxnDefinition(rxnIdx4(i)) = strcat(def_temp(2),{' => '},def_temp(1));
                        end
                    elseif any(idx2)
                        model_coef = model_coef(idx);
                        SEED_coef = SEED_coef(idx2(idx2~=0));
                        if isequal(SEED_coef,model_coef)
                            %no change to rxnDirection
                        elseif isequal(sign(SEED_coef.*(-1)),sign(model_coef)) && strcmp(rxnDirection(rxnIdx4(i)),'>')
                            rxnDirection(rxnIdx4(i)) = {'<'};
                            def_temp = split(rxnDefinition(rxnIdx4(i)), ' => ');
                            rxnDefinition(rxnIdx4(i)) = strcat(def_temp(2),{' <= '},def_temp(1));
                        elseif isequal(sign(SEED_coef.*(-1)),sign(model_coef)) && strcmp(rxnDirection(rxnIdx4(i)),'<')
                            rxnDirection(rxnIdx4(i)) = {'>'};
                            def_temp = split(rxnDefinition(rxnIdx4(i)), ' <= ');
                            rxnDefinition(rxnIdx4(i)) = strcat(def_temp(2),{' => '},def_temp(1));
                        end
                    else
                        warning('%s rxnDirection not checked against %s',string(xref_seed(i,2)),string(xref_seed(i,1)));
                    end
                else
                    [idx,idx2] = ismember(SEED_mets,model_mets);
                    if all(idx2)
                        model_coef = model_coef(idx2);
                        if isequal(SEED_coef,model_coef)
                            %no change to rxnDirection
                        elseif isequal(sign(SEED_coef.*(-1)),sign(model_coef)) && strcmp(rxnDirection(rxnIdx4(i)),'>')
                            rxnDirection(rxnIdx4(i)) = {'<'};
                            def_temp = split(rxnDefinition(rxnIdx4(i)), ' => ');
                            rxnDefinition(rxnIdx4(i)) = strcat(def_temp(2),{' <= '},def_temp(1));
                        elseif isequal(sign(SEED_coef.*(-1)),sign(model_coef)) && strcmp(rxnDirection(rxnIdx4(i)),'<')
                            rxnDirection(rxnIdx4(i)) = {'>'};
                            def_temp = split(rxnDefinition(rxnIdx4(i)), ' <= ');
                            rxnDefinition(rxnIdx4(i)) = strcat(def_temp(2),{' => '},def_temp(1));
                        end
                    elseif any(idx2)
                        SEED_coef = SEED_coef(idx);
                        model_coef = model_coef(idx2(idx2~=0));
                        if isequal(SEED_coef,model_coef)
                            %no change to rxnDirection
                        elseif isequal(sign(SEED_coef.*(-1)),sign(model_coef)) && strcmp(rxnDirection(rxnIdx4(i)),'>')
                            rxnDirection(rxnIdx4(i)) = {'<'};
                            def_temp = split(rxnDefinition(rxnIdx4(i)), ' => ');
                            rxnDefinition(rxnIdx4(i)) = strcat(def_temp(2),{' <= '},def_temp(1));
                        elseif isequal(sign(SEED_coef.*(-1)),sign(model_coef)) && strcmp(rxnDirection(rxnIdx4(i)),'<')
                            rxnDirection(rxnIdx4(i)) = {'>'};
                            def_temp = split(rxnDefinition(rxnIdx4(i)), ' <= ');
                            rxnDefinition(rxnIdx4(i)) = strcat(def_temp(2),{' => '},def_temp(1));
                        end
                    else
                        warning('%s rxnDirection not checked against %s',string(xref_seed(i,2)),string(xref_seed(i,1)));
                    end
                end
            end
            if deltaG > 0
                if ~isequal(seed_rxnInfo{rxnIdx4(i),16},'10000000.0') %check if database contains valid deltaG value
                    mod_deltaG = deltaG - error;
                    if mod_deltaG > 30
                        arrayidx = find(cellfun('isempty',bwd_SEED),1);
                        bwd_SEED(arrayidx,:) = [xref_seed(i,2),deltaG,error,lnRevIdx,rxnDirection(rxnIdx4(i))];
                    else
                        arrayidx = find(cellfun('isempty',rev_SEED),1);
                        rev_SEED(arrayidx,:) = [xref_seed(i,2),deltaG,error,lnRevIdx,rxnDirection(rxnIdx4(i))];
                    end
                else
                    if isequal(rxnDirection(rxnIdx4(i)),'<')
                        arrayidx = find(cellfun('isempty',bwd_SEED),1);
                        bwd_SEED(arrayidx,:) = [xref_seed(i,2),{[]},{[]},{[]},rxnDirection(rxnIdx4(i))];
                    elseif isequal(rxnDirection(rxnIdx4(i)),'=')
                        arrayidx = find(cellfun('isempty',rev_SEED),1);
                        rev_SEED(arrayidx,:) = [xref_seed(i,2),{[]},{[]},{[]},rxnDirection(rxnIdx4(i))];
                    elseif isequal(rxnDirection(rxnIdx4(i)),'>')
                        arrayidx = find(cellfun('isempty',fwd_SEED),1);
                        fwd_SEED(arrayidx,:) = [xref_seed(i,2),{[]},{[]},{[]},rxnDirection(rxnIdx4(i))];
                    else
                        arrayidx = find(cellfun('isempty',no_deltaG3),1);
                        no_deltaG3 = [xref_seed(i,2),{[]},{[]},{[]},{''}];
                    end
                end
            elseif deltaG < 0
                mod_deltaG = deltaG + error;
                if mod_deltaG < -30
                    arrayidx = find(cellfun('isempty',fwd_SEED),1);
                    fwd_SEED(arrayidx,:) = [xref_seed(i,2),deltaG,error,lnRevIdx,rxnDirection(rxnIdx4(i))];
                else
                    arrayidx = find(cellfun('isempty',rev_SEED),1);
                    rev_SEED(arrayidx,:) = [xref_seed(i,2),deltaG,error,lnRevIdx,rxnDirection(rxnIdx4(i))];
                end
            else
                if ~isequal(seed_rxnInfo{rxnIdx4(i),16},'10000000.0') %check if database contains valid deltaG value
                    arrayidx = find(cellfun('isempty',rev_SEED),1);
                    rev_SEED(arrayidx,:) = [xref_seed(i,2),deltaG,error,lnRevIdx,rxnDirection(rxnIdx4(i))];
                else
                    if isequal(rxnDirection(rxnIdx4(i)),'<')
                        arrayidx = find(cellfun('isempty',bwd_SEED),1);
                        bwd_SEED(arrayidx,:) = [xref_seed(i,2),{[]},{[]},{[]},rxnDirection(rxnIdx4(i))];
                    elseif isequal(rxnDirection(rxnIdx4(i)),'=')
                        arrayidx = find(cellfun('isempty',rev_SEED),1);
                        rev_SEED(arrayidx,:) = [xref_seed(i,2),{[]},{[]},{[]},rxnDirection(rxnIdx4(i))];
                    elseif isequal(rxnDirection(rxnIdx4(i)),'>')
                        arrayidx = find(cellfun('isempty',fwd_SEED),1);
                        fwd_SEED(arrayidx,:) = [xref_seed(i,2),{[]},{[]},{[]},rxnDirection(rxnIdx4(i))];
                    end
                end
            end
        else
            arrayidx = find(cellfun('isempty',noMapping3),1);
            noMapping3(arrayidx,:) = [xref_seed(i,2),'-','No MNXRID mapped to rxnSEEDID'];
        end
    end
    
    for i = 1:length(idx_check3)
        arrayidx = find(cellfun('isempty',noMapping3),1);
        noMapping3(arrayidx,:) = [model.rxns(idx_check3(i)),'-','No MNXRID mapped to rxnSEEDID'];
    end
    
    %remove empty cells
    empties = find(cellfun('isempty',bwd_SEED(:,1)));
    bwd_SEED(empties,:) = [];
    empties = find(cellfun('isempty',fwd_SEED(:,1)));
    fwd_SEED(empties,:) = [];
    empties = find(cellfun('isempty',rev_SEED(:,1)));
    rev_SEED(empties,:) = [];
    empties = find(cellfun('isempty',no_deltaG3(:,1)));
    no_deltaG3(empties,:) = [];
    empties = find(cellfun('isempty',noMapping(:,1)));
    noMapping3(empties,:) = [];
    
    %remove downloaded reactions.tsv
    delete('reactions.tsv');
else
    warning('Checking of reaction direction with modelSEED IDs unsuccessful as reactions.tsv is not found');
end

%% Compile checked bwd reactions

bwd_rxns{length(model.rxns),13} = []; %model.rxns grpCont grpCont_uncertainty grpCont_lnRevIdx 
                                      %MetaCyc MetaCyc_rxnDirection MetaCyc_lnRevIdx 
                                      %SEED SEED_uncertainty SEED_lnRevIdx SEED_rxnDirection memote_revIdx original_rxnDirection
if ~isempty(bwd_grpCont)
    rxns = bwd_grpCont(:,1);
    deltaG = bwd_grpCont(:,2);
    uncertainty = bwd_grpCont(:,3);
    lnRevIdx = bwd_grpCont(:,4);
    for i = 1:size(bwd_grpCont,1)
        arrayidx = find(cellfun('isempty',bwd_rxns),1);
        bwd_rxns(arrayidx,1) = rxns(i);
        bwd_rxns(arrayidx,2) = deltaG(i);
        bwd_rxns(arrayidx,3) = uncertainty(i);
        bwd_rxns(arrayidx,4) = lnRevIdx(i);
        [~,idx] = ismember(rxns(i),model.rxns);
        lb = model.lb(idx);
        if lb == 0
            bwd_rxns(arrayidx,13) = {'>'};
        elseif lb == -1000
            bwd_rxns(arrayidx,13) = {'='};
        end
    end
end

if ~isempty(bwd_MetaCyc)
    rxns = bwd_MetaCyc(:,1);
    deltaG = bwd_MetaCyc(:,2);
    lnRevIdx = bwd_MetaCyc(:,3);
    rxnDirection = bwd_MetaCyc(:,4);
    for i = 1:size(bwd_MetaCyc,1)
        if ~isequal(bwd_rxns(1,1),{''})
            temp = find(~cellfun('isempty',bwd_rxns(:,1)));
            temp_rxn = bwd_rxns(temp,1);
            if ismember(rxns(i),temp_rxn)
                [~,arrayidx] =  ismember(rxns(i),temp_rxn);
                bwd_rxns(arrayidx,5) = deltaG(i);
                bwd_rxns(arrayidx,6) = lnRevIdx(i);
                bwd_rxns(arrayidx,7) = rxnDirection(i);
            else
                arrayidx = find(cellfun('isempty',bwd_rxns),1);
                bwd_rxns(arrayidx,1) = rxns(i);
                bwd_rxns(arrayidx,5) = deltaG(i);
                bwd_rxns(arrayidx,6) = lnRevIdx(i);
                bwd_rxns(arrayidx,7) = rxnDirection(i);
                [~,idx] = ismember(rxns(i),model.rxns);
                lb = model.lb(idx);
                if lb == 0
                    bwd_rxns(arrayidx,13) = {'>'};
                elseif lb == -1000
                    bwd_rxns(arrayidx,13) = {'='};
                end
            end
        else
            arrayidx = find(cellfun('isempty',bwd_rxns),1);
            bwd_rxns(arrayidx,1) = rxns(i);
            bwd_rxns(arrayidx,5) = deltaG(i);
            bwd_rxns(arrayidx,6) = lnRevIdx(i);
            bwd_rxns(arrayidx,7) = rxnDirection(i);
            [~,idx] = ismember(rxns(i),model.rxns);
            lb = model.lb(idx);
            if lb == 0
                bwd_rxns(arrayidx,13) =  {'>'};
            elseif lb == -1000
                bwd_rxns(arrayidx,13) = {'='};
            end
        end
    end
end

if ~isempty(bwd_SEED)
    rxns = bwd_SEED(:,1);
    deltaG = bwd_SEED(:,2);
    uncertainty = bwd_SEED(:,3);
    lnRevIdx = bwd_SEED(:,4);
    rxnDirection = bwd_SEED(:,5);
    for i = 1:size(bwd_SEED,1)
        if ~isequal(bwd_rxns(1,1),{''})
            temp = find(~cellfun('isempty',bwd_rxns(:,1)));
            temp_rxn = bwd_rxns(temp,1);
            if ismember(rxns(i),temp_rxn)
                [~,arrayidx] =  ismember(rxns(i),temp_rxn);
                bwd_rxns(arrayidx,8) = deltaG(i);
                bwd_rxns(arrayidx,9) = uncertainty(i);
                bwd_rxns(arrayidx,10) = lnRevIdx(i);
                bwd_rxns(arrayidx,11) = rxnDirection(i);
            else
                arrayidx = find(cellfun('isempty',bwd_rxns),1);
                bwd_rxns(arrayidx,1) = rxns(i);
                bwd_rxns(arrayidx,8) = deltaG(i);
                bwd_rxns(arrayidx,9) = uncertainty(i);
                bwd_rxns(arrayidx,10) = lnRevIdx(i);
                bwd_rxns(arrayidx,11) = rxnDirection(i);
                [~,idx] = ismember(rxns(i),model.rxns);
                lb = model.lb(idx);
                if lb == 0
                    bwd_rxns(arrayidx,13) = {'>'};
                elseif lb == -1000
                    bwd_rxns(arrayidx,13) = {'='};
                end
            end
        else
            arrayidx = find(cellfun('isempty',bwd_rxns),1);
            bwd_rxns(arrayidx,1) = rxns(i);
            bwd_rxns(arrayidx,8) = deltaG(i);
            bwd_rxns(arrayidx,9) = uncertainty(i);
            bwd_rxns(arrayidx,10) = lnRevIdx(i);
            bwd_rxns(arrayidx,11) = rxnDirection(i);
            [~,idx] = ismember(rxns(i),model.rxns);
            lb = model.lb(idx);
            if lb == 0
                bwd_rxns(arrayidx,13) = {'>'};
            elseif lb == -1000
                bwd_rxns(arrayidx,13) = {'='};
            end
        end
    end
end

%remove empty cells
empties = find(cellfun('isempty',bwd_rxns(:,1)));
bwd_rxns(empties,:) = [];

%% Compile checked fwd reactions

fwd_rxns{length(model.rxns),13} = []; %model.rxns grpCont grpCont_uncertainty grpCont_lnRevIdx 
                                      %MetaCyc MetaCyc_rxnDirection MetaCyc_lnRevIdx 
                                      %SEED SEED_uncertainty SEED_lnRevIdx SEED_rxnDirection memote_revIdx original_rxnDirection
if ~isempty(fwd_grpCont)
    rxns = fwd_grpCont(:,1);
    deltaG = fwd_grpCont(:,2);
    uncertainty = fwd_grpCont(:,3);
    lnRevIdx = fwd_grpCont(:,4);
    for i = 1:size(fwd_grpCont,1)
        arrayidx = find(cellfun('isempty',fwd_rxns),1);
        fwd_rxns(arrayidx,1) = rxns(i);
        fwd_rxns(arrayidx,2) = deltaG(i);
        fwd_rxns(arrayidx,3) = uncertainty(i);
        fwd_rxns(arrayidx,4) = lnRevIdx(i);
        [~,idx] = ismember(rxns(i),model.rxns);
        lb = model.lb(idx);
        if lb == 0
            fwd_rxns(arrayidx,13) = {'>'};
        elseif lb == -1000
            fwd_rxns(arrayidx,13) = {'='};
        end
    end
end

if ~isempty(fwd_MetaCyc)
    rxns = fwd_MetaCyc(:,1);
    deltaG = fwd_MetaCyc(:,2);
    lnRevIdx = fwd_MetaCyc(:,3);
    rxnDirection = fwd_MetaCyc(:,4);
    for i = 1:size(fwd_MetaCyc,1)
        if ~isequal(fwd_rxns(1,1),{''})
            temp = find(~cellfun('isempty',fwd_rxns(:,1)));
            temp_rxn = fwd_rxns(temp,1);
            if ismember(rxns(i),temp_rxn)
                [~,arrayidx] =  ismember(rxns(i),temp_rxn);
                fwd_rxns(arrayidx,5) = deltaG(i);
                fwd_rxns(arrayidx,6) = lnRevIdx(i);
                fwd_rxns(arrayidx,7) = rxnDirection(i);
            else
                arrayidx = find(cellfun('isempty',fwd_rxns),1);
                fwd_rxns(arrayidx,1) = rxns(i);
                fwd_rxns(arrayidx,5) = deltaG(i);
                fwd_rxns(arrayidx,6) = lnRevIdx(i);
                fwd_rxns(arrayidx,7) = rxnDirection(i);
                [~,idx] = ismember(rxns(i),model.rxns);
                lb = model.lb(idx);
                if lb == 0
                    fwd_rxns(arrayidx,13) = {'>'};
                elseif lb == -1000
                    fwd_rxns(arrayidx,13) = {'='};
                end
            end
        else
            arrayidx = find(cellfun('isempty',fwd_rxns),1);
            fwd_rxns(arrayidx,1) = rxns(i);
            fwd_rxns(arrayidx,5) = deltaG(i);
            fwd_rxns(arrayidx,6) = lnRevIdx(i);
            fwd_rxns(arrayidx,7) = rxnDirection(i);
            [~,idx] = ismember(rxns(i),model.rxns);
            lb = model.lb(idx);
            if lb == 0
                fwd_rxns(arrayidx,13) =  {'>'};
            elseif lb == -1000
                fwd_rxns(arrayidx,13) = {'='};
            end
        end
    end
end

if ~isempty(fwd_SEED)
    rxns = fwd_SEED(:,1);
    deltaG = fwd_SEED(:,2);
    uncertainty = fwd_SEED(:,3);
    lnRevIdx = fwd_SEED(:,4);
    rxnDirection = fwd_SEED(:,5);
    for i = 1:size(fwd_SEED,1)
        if ~isequal(fwd_rxns(1,1),{''})
            temp = find(~cellfun('isempty',fwd_rxns(:,1)));
            temp_rxn = fwd_rxns(temp,1);
            if ismember(rxns(i),temp_rxn)
                [~,arrayidx] =  ismember(rxns(i),temp_rxn);
                fwd_rxns(arrayidx,8) = deltaG(i);
                fwd_rxns(arrayidx,9) = uncertainty(i);
                fwd_rxns(arrayidx,10) = lnRevIdx(i);
                fwd_rxns(arrayidx,11) = rxnDirection(i);
            else
                arrayidx = find(cellfun('isempty',fwd_rxns),1);
                fwd_rxns(arrayidx,1) = rxns(i);
                fwd_rxns(arrayidx,8) = deltaG(i);
                fwd_rxns(arrayidx,9) = uncertainty(i);
                fwd_rxns(arrayidx,10) = lnRevIdx(i);
                fwd_rxns(arrayidx,11) = rxnDirection(i);
                [~,idx] = ismember(rxns(i),model.rxns);
                lb = model.lb(idx);
                if lb == 0
                    fwd_rxns(arrayidx,13) = {'>'};
                elseif lb == -1000
                    fwd_rxns(arrayidx,13) = {'='};
                end
            end
        else
            arrayidx = find(cellfun('isempty',fwd_rxns),1);
            fwd_rxns(arrayidx,1) = rxns(i);
            fwd_rxns(arrayidx,8) = deltaG(i);
            fwd_rxns(arrayidx,9) = uncertainty(i);
            fwd_rxns(arrayidx,10) = lnRevIdx(i);
            fwd_rxns(arrayidx,11) = rxnDirection(i);
            [~,idx] = ismember(rxns(i),model.rxns);
            lb = model.lb(idx);
            if lb == 0
                fwd_rxns(arrayidx,13) = {'>'};
            elseif lb == -1000
                fwd_rxns(arrayidx,13) = {'='};
            end
        end
    end
end

%remove empty cells
empties = find(cellfun('isempty',fwd_rxns(:,1)));
fwd_rxns(empties,:) = [];

%% Compile checked rev reactions

rev_rxns{length(model.rxns),13} = []; %model.rxns grpCont grpCont_uncertainty grpCont_lnRevIdx 
                                      %MetaCyc MetaCyc_rxnDirection MetaCyc_lnRevIdx 
                                      %SEED SEED_uncertainty SEED_lnRevIdx SEED_rxnDirection memote_revIdx original_rxnDirection
if ~isempty(rev_grpCont)
    rxns = rev_grpCont(:,1);
    deltaG = rev_grpCont(:,2);
    uncertainty = rev_grpCont(:,3);
    lnRevIdx = rev_grpCont(:,4);
    for i = 1:size(rev_grpCont,1)
        arrayidx = find(cellfun('isempty',rev_rxns),1);
        rev_rxns(arrayidx,1) = rxns(i);
        rev_rxns(arrayidx,2) = deltaG(i);
        rev_rxns(arrayidx,3) = uncertainty(i);
        rev_rxns(arrayidx,4) = lnRevIdx(i);
        [~,idx] = ismember(rxns(i),model.rxns);
        lb = model.lb(idx);
        if lb == 0
            rev_rxns(arrayidx,13) = {'>'};
        elseif lb == -1000
            rev_rxns(arrayidx,13) = {'='};
        end
    end
end

if ~isempty(rev_MetaCyc)
    rxns = rev_MetaCyc(:,1);
    deltaG = rev_MetaCyc(:,2);
    lnRevIdx = rev_MetaCyc(:,3);
    rxnDirection = rev_MetaCyc(:,4);
    for i = 1:size(rev_MetaCyc,1)
        if ~isequal(rev_rxns(1,1),{''})
            temp = find(~cellfun('isempty',rev_rxns(:,1)));
            temp_rxn = rev_rxns(temp,1);
            if ismember(rxns(i),temp_rxn)
                [~,arrayidx] =  ismember(rxns(i),temp_rxn);
                rev_rxns(arrayidx,5) = deltaG(i);
                rev_rxns(arrayidx,6) = lnRevIdx(i);
                rev_rxns(arrayidx,7) = rxnDirection(i);
            else
                arrayidx = find(cellfun('isempty',rev_rxns),1);
                rev_rxns(arrayidx,1) = rxns(i);
                rev_rxns(arrayidx,5) = deltaG(i);
                rev_rxns(arrayidx,6) = lnRevIdx(i);
                rev_rxns(arrayidx,7) = rxnDirection(i);
                [~,idx] = ismember(rxns(i),model.rxns);
                lb = model.lb(idx);
                if lb == 0
                    rev_rxns(arrayidx,13) = {'>'};
                elseif lb == -1000
                    rev_rxns(arrayidx,13) = {'='};
                end
            end
        else
            arrayidx = find(cellfun('isempty',rev_rxns),1);
            rev_rxns(arrayidx,1) = rxns(i);
            rev_rxns(arrayidx,5) = deltaG(i);
            rev_rxns(arrayidx,6) = lnRevIdx(i);
            rev_rxns(arrayidx,7) = rxnDirection(i);
            [~,idx] = ismember(rxns(i),model.rxns);
            lb = model.lb(idx);
            if lb == 0
                rev_rxns(arrayidx,13) =  {'>'};
            elseif lb == -1000
                rev_rxns(arrayidx,13) = {'='};
            end
        end
    end
end

if ~isempty(rev_SEED)
    rxns = rev_SEED(:,1);
    deltaG = rev_SEED(:,2);
    uncertainty = rev_SEED(:,3);
    lnRevIdx = rev_SEED(:,4);
    rxnDirection = rev_SEED(:,5);
    for i = 1:size(rev_SEED,1)
        if ~isequal(rev_rxns(1,1),{''})
            temp = find(~cellfun('isempty',rev_rxns(:,1)));
            temp_rxn = rev_rxns(temp,1);
            if ismember(rxns(i),temp_rxn)
                [~,arrayidx] =  ismember(rxns(i),temp_rxn);
                rev_rxns(arrayidx,8) = deltaG(i);
                rev_rxns(arrayidx,9) = uncertainty(i);
                rev_rxns(arrayidx,10) = lnRevIdx(i);
                rev_rxns(arrayidx,11) = rxnDirection(i);
            else
                arrayidx = find(cellfun('isempty',rev_rxns),1);
                rev_rxns(arrayidx,1) = rxns(i);
                rev_rxns(arrayidx,8) = deltaG(i);
                rev_rxns(arrayidx,9) = uncertainty(i);
                rev_rxns(arrayidx,10) = lnRevIdx(i);
                rev_rxns(arrayidx,11) = rxnDirection(i);
                [~,idx] = ismember(rxns(i),model.rxns);
                lb = model.lb(idx);
                if lb == 0
                    rev_rxns(arrayidx,13) = {'>'};
                elseif lb == -1000
                    rev_rxns(arrayidx,13) = {'='};
                end
            end
        else
            arrayidx = find(cellfun('isempty',rev_rxns),1);
            rev_rxns(arrayidx,1) = rxns(i);
            rev_rxns(arrayidx,8) = deltaG(i);
            rev_rxns(arrayidx,9) = uncertainty(i);
            rev_rxns(arrayidx,10) = lnRevIdx(i);
            rev_rxns(arrayidx,11) = rxnDirection(i);
            [~,idx] = ismember(rxns(i),model.rxns);
            lb = model.lb(idx);
            if lb == 0
                rev_rxns(arrayidx,13) = {'>'};
            elseif lb == -1000
                rev_rxns(arrayidx,13) = {'='};
            end
        end
    end
end
cd ../

%remove empty cells
empties = find(cellfun('isempty',rev_rxns(:,1)));
rev_rxns(empties,:) = [];

%Load reversibility index extracted from report in Memote
%Categorize reactions via deltaG from group contribution method
fid2 = fopen('../data/modelCuration/memote_revIdx.tsv');
format = repmat('%s ',1,2);
format = strtrim(format);
rxn_temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',0);
for i = 1:length(rxn_temp)
    rev_index(:,i) = rxn_temp{i};
end
commentLines = startsWith(rev_index(:,1),'#');
rev_index(commentLines,:) = []; % model.rxns memote_revIdx
fclose(fid2);

%add available reversibility index to rev_rxns
[~,idx] = ismember(rev_index(:,1),rev_rxns(:,1));
rev_index(idx==0,:) = [];
idx = idx(idx~=0);

for i = 1:length(idx)
    rev_rxns(idx(i),12) = num2cell(str2double(cell2mat(rev_index(i,2))));
end
%remove from rev_index after adding
[~,idx] = ismember(rev_rxns(:,1),rev_index(:,1));
idx = idx(idx~=0);
rev_index(idx,:) = [];

%add available reversibility index to fwd_rxns
[~,idx] = ismember(rev_index(:,1),fwd_rxns(:,1));
idx = idx(idx~=0);
for i = 1:length(idx)
    fwd_rxns(idx(i),12) = num2cell(str2double(cell2mat(rev_index(i,2))));
end
%remove from rev_index after adding
[~,idx] = ismember(rev_rxns(:,1),rev_index(:,1));
idx = idx(idx~=0);
rev_index(idx,:) = [];

%add available reversibility index to bwd_rxns
[~,idx] = ismember(rev_index(:,1),bwd_rxns(:,1));
idx = idx(idx~=0);
for i = 1:length(idx)
    bwd_rxns(idx(i),12) = num2cell(str2double(cell2mat(rev_index(i,2))));
end
%remove from rev_index after adding
[~,idx] = ismember(rev_rxns(:,1),rev_index(:,1));
idx = idx(idx~=0);
rev_index(idx,:) = [];

%% Further categorization of reactions for manual curation
%Check for potential overlapping reactions e.g. same reaction in model but different databases one with
%deltaG above 30 kJ/mol threshold while the other below 30 kJ/mol threshold

[~,idx] = ismember(bwd_rxns(:,1),rev_rxns(:,1));
conflict_bwd_rev{length(idx(idx~=0)),13} = [];

for i = 1:length(idx)
    if idx(i) ~= 0
        arrayidx = find(cellfun('isempty',conflict_bwd_rev),1);
        conflict_bwd_rev(arrayidx,1) = bwd_rxns(i,1);
        conflict_bwd_rev(arrayidx,13) = bwd_rxns(i,13);
        for j = 2:12
            if ~isequal(bwd_rxns(i,j),{''})
                conflict_bwd_rev(arrayidx,j) = bwd_rxns(i,j);
            end
        end
        for j = 2:12
            if ~isequal(rev_rxns(idx(i),j),{''})
                conflict_bwd_rev(arrayidx,j) = rev_rxns(idx(i),j);
            end
        end
    end
end

%remove reactions with conflict from bwd_rxns and rev_rxns
[~,idx2] = ismember(rev_rxns(:,1),bwd_rxns(:,1));
idx = idx(idx~=0);
rev_rxns(idx,:) = [];
idx2 = idx2(idx2~=0);
bwd_rxns(idx2,:) = [];

%Check for potential overlapping reactions e.g. same reaction in model but different databases one with
%deltaG above 30 kJ/mol threshold while the other below 30 kJ/mol threshold
[~,idx] = ismember(fwd_rxns(:,1),rev_rxns(:,1));
conflict_fwd_rev{length(idx(idx~=0)),13} = [];

for i = 1:length(idx)
    if idx(i) ~= 0
        arrayidx = find(cellfun('isempty',conflict_fwd_rev),1);
        conflict_fwd_rev(arrayidx,1) = fwd_rxns(i,1);
        conflict_fwd_rev(arrayidx,13) = fwd_rxns(i,13);
        for j = 2:12
            if ~isequal(fwd_rxns(i,j),{''})
                conflict_fwd_rev(arrayidx,j) = fwd_rxns(i,j);
            end
        end
        for j = 2:12
            if ~isequal(rev_rxns(idx(i),j),{''})
                conflict_fwd_rev(arrayidx,j) = rev_rxns(idx(i),j);
            end
        end
    end
end

%remove reactions with conflict from fwd_rxns and rev_rxns
[~,idx2] = ismember(rev_rxns(:,1),fwd_rxns(:,1));
idx = idx(idx~=0);
rev_rxns(idx,:) = [];
idx2 = idx2(idx2~=0);
fwd_rxns(idx2,:) = [];

checked_rxns = unique([fwd_rxns(:,1);bwd_rxns(:,1);rev_rxns(:,1)]);

%Categorise reactions based on deltaG values
%Categorise reactions in bwd_rxns
bwdfwd_unmatched{length(bwd_rxns),13} = [];
bwdrev_unmatched{length(bwd_rxns),13} = [];

for i = 1:size(bwd_rxns,1)
    if strcmp(bwd_rxns(i,13),'>')
        arrayidx = find(cellfun('isempty',bwdfwd_unmatched),1);
        bwdfwd_unmatched(arrayidx,1) = bwd_rxns(i,1);
        bwdfwd_unmatched(arrayidx,13) = bwd_rxns(i,13);
        for j = 2:12
            if ~isequal(bwd_rxns(i,j),{''})
                bwdfwd_unmatched(arrayidx,j) = bwd_rxns(i,j);
            end
        end
    elseif strcmp(bwd_rxns(i,13),'=')
        arrayidx = find(cellfun('isempty',bwdrev_unmatched),1);
        bwdrev_unmatched(arrayidx,1) = bwd_rxns(i,1);
        bwdrev_unmatched(arrayidx,13) = bwd_rxns(i,13);
        for j = 2:12
            if ~isequal(bwd_rxns(i,j),{''})
                bwdrev_unmatched(arrayidx,j) = bwd_rxns(i,j);
            end
        end
    end
end

%remove empty cells
empties = find(cellfun('isempty',bwdfwd_unmatched(:,1)));
bwdfwd_unmatched(empties,:) = [];
empties = find(cellfun('isempty',bwdrev_unmatched(:,1)));
bwdrev_unmatched(empties,:) = [];

%Categorise reactions in fwd_rxns
fwd_matched{length(fwd_rxns),13} = [];
fwd_unmatched{length(fwd_rxns),13} = [];

for i = 1:size(fwd_rxns,1)
    if strcmp(fwd_rxns(i,13),'>')
        arrayidx = find(cellfun('isempty',fwd_matched),1);
        fwd_matched(arrayidx,1) = fwd_rxns(i,1);
        fwd_matched(arrayidx,13) = fwd_rxns(i,13);
        for j = 2:12
            if ~isequal(fwd_rxns(i,j),{''})
                fwd_matched(arrayidx,j) = fwd_rxns(i,j);
            end
        end
    elseif strcmp(fwd_rxns(i,13),'=')
        arrayidx = find(cellfun('isempty',fwd_unmatched),1);
        fwd_unmatched(arrayidx,1) = fwd_rxns(i,1);
        fwd_unmatched(arrayidx,13) = fwd_rxns(i,13);
        for j = 2:12
            if ~isequal(fwd_rxns(i,j),{''})
                fwd_unmatched(arrayidx,j) = fwd_rxns(i,j);
            end
        end
    end
end

%remove empty cells
empties = find(cellfun('isempty',fwd_matched(:,1)));
fwd_matched(empties,:) = [];
empties = find(cellfun('isempty',fwd_unmatched(:,1)));
fwd_unmatched(empties,:) = [];

%Categorise reactions in rev_rxns
rev_matched{length(rev_rxns),13} = [];
rev_unmatched{length(rev_rxns),13} = [];

for i = 1:size(rev_rxns,1)
    if strcmp(rev_rxns(i,13),'>')
        arrayidx = find(cellfun('isempty',rev_unmatched),1);
        rev_unmatched(arrayidx,1) = rev_rxns(i,1);
        rev_unmatched(arrayidx,13) = rev_rxns(i,13);
        for j = 2:12
            if ~isequal(rev_rxns(i,j),{''})
                rev_unmatched(arrayidx,j) = rev_rxns(i,j);
            end
        end
    elseif strcmp(rev_rxns(i,13),'=')
        arrayidx = find(cellfun('isempty',rev_matched),1);
        rev_matched(arrayidx,1) = rev_rxns(i,1);
        rev_matched(arrayidx,13) = rev_rxns(i,13);
        for j = 2:12
            if ~isequal(rev_rxns(i,j),{''})
                rev_matched(arrayidx,j) = rev_rxns(i,j);
            end
        end
    end
end

%remove empty cells
empties = find(cellfun('isempty',rev_matched(:,1)));
rev_matched(empties,:) = [];
empties = find(cellfun('isempty',rev_unmatched(:,1)));
rev_unmatched(empties,:) = [];

%Complie all rxns in one variable
allrxns = [rev_matched;rev_unmatched;fwd_matched;fwd_unmatched;...
    bwdfwd_unmatched;bwdrev_unmatched;conflict_fwd_rev;conflict_bwd_rev];

clearvars -except allrxns

%{
    This segment was used to generate rxnDirectionInfo.tsv, which was
    subsequently edited and used in changerxnDirection.m:
    
    file = allrxns';
    fid = fopen('rxnDirectionInfo.tsv', 'wt');
    format = repmat('%s\t',1,13);
    format(end) = 'n';
    format2 = ['%s\t',repmat('%5.4f\t',1,5),'%s\t',repmat('%5.4f\t',1,3),'%s\t','%5.4f\t','%s\n'];
    fprintf(fid, format,'#model.rxns','grpCont_deltaG','grpCont_uncertainty','grpCont_lnRevIdx',...
        'MetaCyc_deltaG','MetaCyc_rxnDirection','MetaCyc_lnRevIdx',...
        'SEED_deltaG','SEED_uncertainty','SEED_lnRevIdx','SEED_rxnDirection','memote_revIdx','model_rxnDirection');
    fprintf(fid, format2,file{:});
    fclose(fid);
%}

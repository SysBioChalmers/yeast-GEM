%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SubstrateUsage
%
% Automatically adds exchange reactions for every metabolite in 
% data/physiology/Biolog_substrate.tsv and checks whether
% it can be used as a "solo" substrate.
%
% NOTE: requires COBRA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load model:
cd ..
model = loadYeastModel;

% Load data:
fid2 = fopen('../data/physiology/Biolog_substrate.tsv');
substrate = textscan(fid2,'%s %s %s %s %s %s','Delimiter','\t','HeaderLines',1);
SubBiologName = substrate{1};
SubModelName  = substrate{2};
Subtype       = substrate{3};
BiologUsage   = substrate{4};
fclose(fid2);

% Main loop:
cd otherChanges
ExchRxn             = '';
TransRxn            = '';
GapfillMets         = '';
TransEotherfillMets = '';
TransECfillMets     = '';
FBAresult           = '';
for i = 1:length(SubBiologName)
    newModel = model;
    if ~isequal(SubModelName{i},'')
        metE = strcat(SubModelName{i},' [extracellular]');
        metC = strcat(SubModelName{i},' [cytoplasm]');
        [~,metEOrd] = ismember(metE,newModel.metNames);
        [~,metCOrd] = ismember(metC,newModel.metNames);
        %check whther exchange reaction for this met exists
        duplicateE = find(sum(newModel.S ~= 0, 1) == 1 & any(newModel.S == -1, 1) & any(newModel.S(metEOrd(metEOrd~=0), :), 1));
        if ~isempty(duplicateE)
            warning(['Model already has the same exchange reaction you tried to add: ', newModel.rxns{duplicateE}]);
            ExchRxn  = newModel.rxns{duplicateE};
            ExchLB   = newModel.lb(duplicateE);
        else
            if metEOrd == 0
                newID       = getNewIndex(newModel.mets);
                metE_ID     = strcat('s_',newID,'[e]');
                newModel    = addMetabolite(newModel,metE_ID,'metName',metE);
                [~,metEOrd] = ismember(metE,newModel.metNames);
            end
            %adding exchange rxn for this metE
            newID    = getNewIndex(newModel.rxns);
            ExchRxn  = ['r_' newID];
            ExchLB   = 0;
        end
        metE_ID = newModel.mets{metEOrd};
        [newModel,rxnIDexists] = addReaction(newModel,ExchRxn, ...
            'reactionName', [SubModelName{i}, ' exchange'], ...
            'metaboliteList', cellstr(metE_ID), 'stoichCoeffList', -1, ...
            'lowerBound', ExchLB, 'upperBound', 1000, 'subSystem', '', ...
            'checkDuplicate', false);
        newModel = rmfield(newModel,'grRules');
        
        % Fix confidence score:
        SubRxnIndex = findRxnIDs(newModel,ExchRxn);
        newModel.rxnConfidenceScores(SubRxnIndex) = NaN;    %exchange rxns
        newModel.rxnNotes{SubRxnIndex} = ['NOTES: added after the Biolog update (PR #149); ', newModel.rxnNames{SubRxnIndex}];
        
        % Change media:
        newModel_test = newModel;
        exchangeRxns  = findExcRxns(newModel_test);
        newModel_test.lb(exchangeRxns) = 0;
        newModel_test.ub(exchangeRxns) = 1000;
        commonExchanges = {...
            'r_1992'; ... % oxygen exchange
            'r_1861'; ... % iron exchange
            'r_1832'; ... % hydrogen exchange
            };
        glcExchange = {'r_1714'}; % D-glucose exchange
        amoExchange = {'r_1654'}; % ammonium exchange
        phoExchange = {'r_2005'}; % phosphate exchange
        sulExchange = {'r_2060'}; % phosphate exchange
        uptakeRxnIndexes = findRxnIDs(newModel_test,commonExchanges);
        amoExchangeIndex = findRxnIDs(newModel_test,amoExchange);
        phoExchangeIndex = findRxnIDs(newModel_test,phoExchange);
        glcExchangeIndex = findRxnIDs(newModel_test,glcExchange);
        sulExchangeIndex = findRxnIDs(newModel_test,sulExchange);
        if length(uptakeRxnIndexes) ~= 3
            error('Not all exchange reactions were found.')
        end
        newModel_test.lb(uptakeRxnIndexes) = -1000;
        newModel_test.lb(glcExchangeIndex) = -10;
        newModel_test.lb(amoExchangeIndex) = -1000;
        newModel_test.lb(phoExchangeIndex) = -1000;
        newModel_test.lb(sulExchangeIndex) = -1000;
        if Subtype{i} == 'C'
            newModel_test.lb(glcExchangeIndex) = 0;
        elseif Subtype{i} == 'N'
            newModel_test.lb(amoExchangeIndex) = 0;
        elseif Subtype{i} == 'P'
            newModel_test.lb(phoExchangeIndex) = 0;
        elseif Subtype{i} == 'S'
            newModel_test.lb(sulExchangeIndex) = 0;
        end
        newModel_test.lb(SubRxnIndex) = -10;
        
        % Simulate model:
        sol = optimizeCbModel(newModel_test);
        if sol.obj > 0.000001
            model = newModel;
            fprintf(metE,'can be used as solo substrate\n');
            FBAresult = [FBAresult;SubBiologName(i,1),SubModelName(i,1),Subtype(i,1),BiologUsage(i,1),'G'];
        else
            FBAresult = [FBAresult;SubBiologName(i,1),SubModelName(i,1),Subtype(i,1),BiologUsage(i,1),'NG'];
            if strcmp(BiologUsage(i,1),'G')
                fprintf(metE,'cannot be used as solo substrate, pathways need to be manually checked\n');
                if numel(find(newModel_test.S(metEOrd,:))) < 2 && metCOrd~=0 %metE appear only once without metC
                    warning(['Transport reaction should be added for:', newModel_test.metNames{metEOrd}])
                    TransECfillMets = [TransECfillMets;newModel_test.metNames(metEOrd)];
                elseif numel(find(newModel_test.S(metEOrd,:))) < 2 && metCOrd ==0 %metE appear only once without metC
                    S = regexp(newModel_test.metNames{metEOrd}, ' [', 'split');
                    met_temp = char(S(1));
                    metEMapping = find(strncmp(met_temp,model.metNames,length(met_temp)));
                    if isempty(metEMapping)
                        warning(['mets only appear once, fill the gap for mets:', newModel_test.metNames{metEOrd}])
                        GapfillMets = [GapfillMets;newModel_test.metNames(metEOrd)];
                    else
                        TransEotherfillMets = [TransEotherfillMets;newModel_test.metNames(metEOrd);model.metNames(metEMapping)];
                    end
                end
            end
        end
    else
        FBAresult = [FBAresult;SubBiologName(i,1),SubModelName(i,1),Subtype(i,1),BiologUsage(i,1),'NG'];
    end
end

% Update results from model:
cd ..
fid2 = fopen('../data/physiology/Biolog_substrate.tsv','w');
formatSpec = '%s\t%s\t%s\t%s\t%s\n';
fprintf(fid2,formatSpec,'Substrate','Name_in_Model','Substrate_type','Growth_Biolog','Growth_Model');
for i = 1:length(FBAresult)
    fprintf(fid2,formatSpec,char(FBAresult(i,1)),char(FBAresult(i,2)),char(FBAresult(i,3)),char(FBAresult(i,4)),char(FBAresult(i,5)));
end

% Save model:
saveYeastModel(model)
cd modelCuration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

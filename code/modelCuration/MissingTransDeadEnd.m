function [model,canbesolved,DEM,DEM_reason] = MissingTransDeadEnd(model)

% This function is to detect whether the deadend metabolites can be solved by adding a transport reaction.
% Output: canbesolved is a list of deadend metabolites that can be solved by adding a transport rxn.
%         DEM is a list of deadend metaolites 
%         DEM_reason is a combined list with possible reasons for all
%         deadend mets.
%         output format for canbesolved: deadentMetsName MNXID deadendMetsInanotherComps Possible solution

model_org = model;
changeCobraSolver('gurobi', 'LP');
exchangeRxns = findExcRxns(model);
model.lb(exchangeRxns) = -1000;
model.ub(exchangeRxns) = 1000;
% DEM are list of total deadendmets in the model
DEM_Idx = detectDeadEnds(model);
DEM = model.metNames(DEM_Idx);
DEMMNX = model.metMetaNetXID(DEM_Idx);

model_r = ravenCobraWrapper(model);

DEM_nosolve = {'deadentMetsName','MNXID','deadendMetsInanotherComps','Possible solution'};
% DEM_nosolve is a list contains all deadend mets that can not be solved by adding a transport in the model
DEM_maysolve = {'deadentMetsName','MNXID','deadendMetsInanotherComps','Possible solution'};
% DEM_maysolve refers to deadend mets that may be solved by adding a transport
% rxn or change lb, which will be further checked in the later part of this
% code.
for i = 1:length(DEM)
    mets_Idx = find(strcmp(model_r.metNames,model_r.metNames(DEM_Idx(i))));
    metsinOtherComps_Idx = setdiff(mets_Idx,DEM_Idx(i));
    metinOtherExE_Idx = metsinOtherComps_Idx(model_r.metComps(metsinOtherComps_Idx) ~= 3);% 3 refers to extracelluar compartment in model_r.metComps
    if ~isempty(metinOtherExE_Idx)
        for j = 1:length(metinOtherExE_Idx)
            if isempty(find(ismember(metinOtherExE_Idx(j),DEM_Idx), 1))
                transrxn = intersect(find(model.S(DEM_Idx,:)~=0),find(model.S(metinOtherExE_Idx(j),:)~=0));
                if ~isempty(transrxn)
                    for m = 1:length(transrxn)
                        if model.lb(transrxn(m)) == 0
                            DEM_maysolve = [DEM_maysolve;DEM(i),DEMMNX(i),model.metNames(metinOtherExE_Idx(j)), ['change lb for rxn ',model.rxns{transrxn(m)}]];
                        elseif model.lb(transrxn) == -1000
                            DEM_nosolve = [DEM_nosolve;DEM(i),DEMMNX(i),model.metNames(metinOtherExE_Idx(j)), ['has transport reaction ', model.rxns{transrxn(m)},' for ', model.metNames{metinOtherExE_Idx(j)}, ' but not work']];
                        end
                    end
                else
                    if ~isempty(find(model_r.metComps(metinOtherExE_Idx(j)) == 1))
                        DEM_maysolve = [DEM_maysolve;DEM(i),DEMMNX(i),model.metNames(metinOtherExE_Idx(j)), ['add a transport reaction from cytosol for ',model.metNames{metinOtherExE_Idx(j)}]];
                    else
                        DEM_maysolve = [DEM_maysolve;DEM(i),DEMMNX(i),model.metNames(metinOtherExE_Idx(j)), ['add a transport reaction for ',model.metNames{metinOtherExE_Idx(j)}]];
                    end
                end
            else
                DEM_nosolve = [DEM_nosolve;DEM(i),DEMMNX(i),model.metNames(metinOtherExE_Idx(j)),' mets in other compartment is also deadend'];
            end
        end
    else
        DEM_nosolve = [DEM_nosolve;DEM(i),DEMMNX(i),'noMetsInOtherComps','mets only appear in one compartment'];
    end
end
DEM_reason = [DEM_maysolve; DEM_nosolve];

% trying to add a transport reaction to see whether the deadend metaolites
% canbesolved is a list that deadend mets can be linked into the model by adding a transport rxns
canbesolved = DEM_maysolve(1,:); % title line
for i = 2:length(DEM_maysolve(:,1))% There is a title line
    if ~isempty(cell2mat(DEM_maysolve(i,1)))
        if strncmpi(DEM_maysolve(i,4),'add a transport reaction from cytosol',37)
            mets = [DEM_maysolve(i,1),DEM_maysolve(i,3)];
            [~,metindex] = ismember(mets,model.metNames);
            metsID = model.mets(metindex);
            cd ../otherChanges
            newID    = getNewIndex(model.rxns);
            cd ../modelCuration/
            TransRxn  = ['r_' newID];
            newModel = addReaction(model,TransRxn, ...
                'reactionName', [DEM_maysolve{i,2}, ' transport'], ...
                'metaboliteList', metsID, 'stoichCoeffList', [-1 1], ...
                'lowerBound', -1000, 'upperBound', 1000, 'subSystem', '', ...
                'checkDuplicate', false);
            DEM_Idx_temp = detectDeadEnds(newModel);
            if ~ismember(DEM_maysolve{i,1},model.metNames(DEM_Idx_temp)) % deadmet is not in the list
                canbesolved = [canbesolved;DEM_maysolve(i,:)];
                model = newModel;
            end
        elseif strncmpi(DEM_maysolve(i,4),'change lb',9)
            rxn_temp = DEM_maysolve{i,4};
            rxnID = rxn_temp(end-5:end);
            [~,rxnindex] = ismember({rxnID},model.rxns);
            newModel = model;
            newModel.lb(rxnindex) = -1000;
            DEM_Idx_temp = detectDeadEnds(newModel);
           if ~ismember(DEM_maysolve{i,1},model.metNames(DEM_Idx_temp))
                canbesolved = [canbesolved;DEM_maysolve(i,:)];
                model = newModel;
           end
        end
    end
end

% return the origninal model
model = model_org;

% clear redundant variables 
clearvars -except DEM DEM_maysolve DEM_nosolve DEM_reason canbesolved model
end

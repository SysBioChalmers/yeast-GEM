% This function is to identify charge and mass balance for reactions
%
% Input: model and reaction identifier(s) in model.rxns (e.g. 'r_0001')
% Output: MassChargeresults
%         For exchange reactions:
%         reaction identifier(s) and 'exchange' will be displayed in
%         MassChargeresults
%
%         For balanced reactions:
%         reaction identifier(s) and 'pass' will be displayed in
%         MassChargeresults
%
%         For reaction identifier(s) which does not match model.rxns: 
%         'rxnID not found' and 'manual check required' will be displayed
%         in MassChargeresults
%
%         For reactions with elemental and/or charge imbalance(s):  
%         reaction identifier(s), specific imbalance issue(s) to check, 
%         individual charges of each metabolite in the reaction (if there 
%         is charge difference), charge difference and elemental difference
%         will be displayed in MassChargeresults
%
% NOTE: getElementalBalance.m is a function from RAVEN
%
% Usage: MassChargeresults = checkRxnBalance(model,rxn)
%
% modified from 'checkBalanceforSce.m'

function MassChargeresults = checkRxnBalance(model,rxn)
exchangeRxns = findExcRxns(model);
MassChargeresults{length(rxn),5} = [];
for i = 1:length(rxn)
    [~,rxnID] = ismember(rxn(i),model.rxns);
    if rxnID ~= 0
        if exchangeRxns(rxnID) == 1
            arrayidx = find(cellfun('isempty', MassChargeresults),1);
            MassChargeresults(arrayidx,1) = model.rxns(rxnID);
            MassChargeresults(arrayidx,2) = {'exchange'};
        else
            %check mass balance
            balanceStructure=getElementalBalance(model,rxn(i));
            unbalancedidx = find(balanceStructure.balanceStatus==0);
            missingidx = find(balanceStructure.balanceStatus==-1);
            idx = union(unbalancedidx,missingidx);
            for j=1:length(idx)
                dif = balanceStructure.leftComp(idx(j),:) - balanceStructure.rightComp(idx(j),:);
                mets = find(~dif==0);
                difference = '';
                for k=1:length(mets)
                    difference = strcat(difference,balanceStructure.elements.abbrevs(mets(k)));
                    difference = strcat(difference,num2str(dif(mets(k))));
                end
                out = difference;
            end
            %check charge balance
            mets=find(model.S(:,rxnID));
            coef = model.metCharges(mets);
            if length(mets) == length(coef)
                indvCharge = model.S(mets,rxnID).*coef;
                balanceCharge = sum(model.S(mets,rxnID).*coef);
            else
                balanceCharge = -1;
            end
            if balanceStructure.balanceStatus == 1 && balanceCharge == 0
                arrayidx = find(cellfun('isempty', MassChargeresults),1);
                MassChargeresults(arrayidx,1) = model.rxns(rxnID);
                MassChargeresults(arrayidx,2) = {'pass'};
            elseif balanceStructure.balanceStatus == 1 && balanceCharge ~= 0
                arrayidx = find(cellfun('isempty', MassChargeresults),1);
                MassChargeresults(arrayidx,1) = model.rxns(rxnID);
                MassChargeresults(arrayidx,2) = {'checkMetCharge'};
                MassChargeresults(arrayidx,3) = {indvCharge};
                MassChargeresults(arrayidx,4) = {balanceCharge};
            elseif balanceStructure.balanceStatus == 0 && balanceCharge == 0
                arrayidx = find(cellfun('isempty', MassChargeresults),1);
                MassChargeresults(arrayidx,1) = model.rxns(rxnID);
                MassChargeresults(arrayidx,2) = {'checkMetFormula'};
                MassChargeresults(arrayidx,5) = out;
            elseif balanceStructure.balanceStatus == 0 && balanceCharge ~= 0
                arrayidx = find(cellfun('isempty', MassChargeresults),1);
                MassChargeresults(arrayidx,1) = model.rxns(rxnID);
                MassChargeresults(arrayidx,2) = {'checkMetChargeAndFormula'};
                MassChargeresults(arrayidx,3) = {indvCharge};
                MassChargeresults(arrayidx,4) = {balanceCharge};
                MassChargeresults(arrayidx,5) = out;
            elseif balanceStructure.balanceStatus == -1
                arrayidx = find(cellfun('isempty', MassChargeresults),1);
                MassChargeresults(arrayidx,1) = model.rxns(rxnID);
                MassChargeresults(arrayidx,2) = {'unbalanced due to missing information'};
                MassChargeresults(arrayidx,3) = {indvCharge};
                MassChargeresults(arrayidx,4) = {balanceCharge};
                MassChargeresults(arrayidx,5) = {'-'};
            end
        end
    else
        arrayidx = find(cellfun('isempty', MassChargeresults),1);
        MassChargeresults(arrayidx,1) = {'rxnID not found'};
        MassChargeresults(arrayidx,2) = {'manual check required'};
    end
end

end

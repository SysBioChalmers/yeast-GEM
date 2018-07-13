%This function is to identify charge and mass balance for new reactions
%the input is model and new reaction ID lists
%output is a result with pass or error


% NOTE: getElementalBalance.m is a function from raven
%
% Feiran Li
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MassChargeresults] = CheckBalanceforSce(model,rxn,MassChargeresults) 

exchangeRxns = findExcRxns(model);
%MassChargeresults = {};
for i = 1:length(rxn)
[~,rxnID] = ismember(rxn(i),model.rxns);
met = (find(model.S(:,rxnID)) ~= 0);
%[nRxn,nMet] = size(rxn,met);
if exchangeRxns(rxnID) == 1
  MassChargeresults = [MassChargeresults; model.rxns(rxnID),'exchange'];
else  
    %check mass balance
    balanceStructure=getElementalBalance(model,rxn(i));
    %check charge balance
    mets=find(model.S(:,rxnID));
    coef = model.metCharges(mets);
    if length(mets) == length(coef)
    balanceCharge = sum(model.S(mets,rxnID).*coef);
    else
        balanceCharge = -1;
    end
    if balanceStructure.balanceStatus == 1 && balanceCharge == 0 
        MassChargeresults = [MassChargeresults; model.rxns(rxnID),'pass'];
    elseif balanceStructure.balanceStatus == 1 && balanceCharge == -1
        MassChargeresults = [MassChargeresults; model.rxns(rxnID),'checkMetCharge'];
    elseif balanceStructure.balanceStatus ~= -1 && balanceCharge == 0
        MassChargeresults = [MassChargeresults; model.rxns(rxnID),'checkMetFormula'];  
    elseif balanceStructure.balanceStatus ~= -1 && balanceCharge == -1
        MassChargeresults = [MassChargeresults; model.rxns(rxnID),'checkMetChargeAndFormula']; 
    else
         MassChargeresults = [MassChargeresults; model.rxns(rxnID),'error'];
end
end
end

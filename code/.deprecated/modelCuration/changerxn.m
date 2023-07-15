function model = changerxn(model,rxnID,rxnformula,grRule)
% changerxn
%   Change the metabolites in the reaction
%
%   model           COBRA model structure
%   rxnID           ID of the reaction i.e. identifier in model.rxns
%   rxnFormula      Reaction formula in string format 
%                   (A [comp] + B [comp] -> C [comp])
%                   For metabolites with space in model.metNames e.g. 
%                   1-pyrroline-3-hydroxy-5-carboxylic acid, use symbol '&'
%                   to substitute space ' ' i.e. 1-pyrroline-3-hydroxy-5-carboxylic&acid
%   grRule          Gene-reaction rule in boolean format (and/or allowed)
%                   (opt, default model.grRule from RAVEN model structure)               
%
%   Note: model.mets is generated for new metabolites with prefix 's_'
%
%   Usage: model = changerxn(model,rxnID,rxnformula,grRule)
%

[~,idx] = ismember(rxnID,model.rxns); 
if nargin < 4
    if isfield(model,'grRules')
        grRule = model.grRules{idx};
    else
        model1 = ravenCobraWrapper(model);
         grRule = model1.grRules{idx};       
    end
end

rxnformula = strrep(rxnformula,' [','[');
[metaboliteList, stoichCoeffList, revFlag] = parseRxnFormula(rxnformula);
metaboliteList = strrep(metaboliteList,'[',' [');
metaboliteList = strrep(metaboliteList,'&',' ');
comps = split(metaboliteList', ' [');
comps = comps(:,2);
comps = strrep(comps,']','');
CONValldata = cat(2,model.compNames,model.comps);
[~,b] = ismember(comps,CONValldata(:,1));
comps = CONValldata(b,2);

%mapping mets to model.metnames, get s_ index for new mets
for j = 1:length(metaboliteList)
    [~,metindex] = ismember(metaboliteList(j),model.metNames);
    if metindex ~= 0
        mets(j) = model.mets(metindex);
    elseif metindex == 0
        cd ../otherChanges/
        newID = getNewIndex(model.mets);
        cd ../modelCuration/
        mets(j) = strcat('s_',newID,'[',comps(j),']');
        model = addMetabolite(model,char(mets(j)), ...
                            'metName',metaboliteList{j});
    end
end

    
[model, rxnIDexists] = addReaction(model,...
                                    rxnID,...
                                    'reactionName', model.rxnNames{idx},...
                                    'metaboliteList',mets,...
                                    'stoichCoeffList',stoichCoeffList,...
                                    'reversible',revFlag,...
                                    'geneRule',grRule,...
                                    'checkDuplicate',1);
 
end

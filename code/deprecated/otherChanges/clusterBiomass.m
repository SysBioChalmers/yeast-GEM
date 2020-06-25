%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = clusterBiomass(model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = clusterBiomass(model)

%Introduce 4 new pseudo-metabolites: protein, carbohydrate, RNA and DNA
model   = addPseudoMet(model,'protein');
model   = addPseudoMet(model,'carbohydrate');
model   = addPseudoMet(model,'RNA');
model   = addPseudoMet(model,'DNA');
protRxn = double(strcmp(model.metNames,'protein [cytoplasm]'));
carbRxn = double(strcmp(model.metNames,'carbohydrate [cytoplasm]'));
rnaRxn  = double(strcmp(model.metNames,'RNA [cytoplasm]'));
dnaRxn  = double(strcmp(model.metNames,'DNA [cytoplasm]'));

%The original biomass pseudo-rxn starts with only the new pseudo-mets as substrates:
bioRxn = -(protRxn + carbRxn + rnaRxn + dnaRxn);

%Go through biomass pseudo-rxn and re-assign each element to where it should:
bioPos = strcmp(model.rxns,'r_4041');
for i = 1:length(model.mets)
    Six = model.S(i,bioPos);
    if Six ~= 0        
        name   = model.metNames{i};
        isProt = contains(name,'tRNA');
        isCarb = sum(strcmpi({'(1->3)-beta-D-glucan [cell envelope]', ...
                              '(1->6)-beta-D-glucan [cell envelope]', ...
                              'chitin [cytoplasm]','glycogen [cytoplasm]', ...
                              'mannan [cytoplasm]','trehalose [cytoplasm]'},name)) == 1;
        isRNA = sum(strcmpi({'AMP [cytoplasm]','CMP [cytoplasm]', ...
                             'GMP [cytoplasm]','UMP [cytoplasm]'},name)) == 1;
        isDNA = sum(strcmpi({'dAMP [cytoplasm]','dCMP [cytoplasm]', ...
                             'dGMP [cytoplasm]','dTMP [cytoplasm]'},name)) == 1;
        if isProt
            protRxn(i) = Six;
        elseif isCarb
            carbRxn(i) = Six;
        elseif isRNA
            rnaRxn(i) = Six;
        elseif isDNA
            dnaRxn(i) = Six;
        else
            bioRxn(i) = Six;
        end
    end
end

%Add the new reactions to the model:
model = addPseudoRxn(model,[],'protein pseudoreaction',protRxn);
model = addPseudoRxn(model,[],'carbohydrate pseudoreaction',carbRxn);
model = addPseudoRxn(model,[],'RNA pseudoreaction',rnaRxn);
model = addPseudoRxn(model,[],'DNA pseudoreaction',dnaRxn);
model = addPseudoRxn(model,'r_4041','biomass pseudoreaction',bioRxn);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = addPseudoMet(model,metName)

metID  = ['s_' getNewIndex(model.mets) '[c]'];
model  = addMetabolite(model,metID,'metName',[metName ' [cytoplasm]']);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = addPseudoRxn(model,rxnID,rxnName,stoich)

if isempty(rxnID)
    rxnID = ['r_' getNewIndex(model.rxns)];
end

model = addReaction(model, ...                      %model
                   {rxnID,rxnName}, ...             %rxn
                    model.mets(stoich ~= 0), ...	%metabolites
                    stoich(stoich ~= 0),     ...    %stoichiometry
                    false, ...                      %reversibility
                    0, ...                          %LB
                    1000, ...                       %UB
                    0);                             %c

printRxnFormula(model,rxnID,true,true,true);
                
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
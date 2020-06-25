%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = addConfidenceScores(model)
% Rough confidence scores for reaction
% Reactions with pubmedID and with gene information: 3
% Reactions with gene but without pubmedID: 2 
% Reactions without gene but need for modelling: 1
% Reactions without gene: 0
% Exchange reactions: NaN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = addConfidenceScores(model)

rxnConfidenceScores = zeros(size(model.rxns));

for i = 1:length(model.rxns)
    if ~isempty(model.rules{i})
        if ~isempty(model.rxnReferences{i})
            rxnConfidenceScores(i) = 3;
        else
            rxnConfidenceScores(i) = 2;
        end
    else
        rxnName = model.rxnNames{i};
        rxnNotes = model.rxnNotes{i};
        if contains(rxnName,'exchange')
            rxnConfidenceScores(i) = NaN;
        elseif contains(rxnName,'SLIME rxn') || contains(rxnName,'pseudoreaction')  || contains(rxnNotes,'Biolog update') || contains(rxnNotes,'BiomassUpdate')
            rxnConfidenceScores(i) = 1;
        else
            metNames = model.metNames(model.S(:,i) ~= 0);
            for j = 1:length(metNames)
                if contains(metNames{j},'backbone [') || contains(metNames{j},'chain [')
                    rxnConfidenceScores(i) = 1;
                end
            end
        end
    end
end

model.rxnConfidenceScores = rxnConfidenceScores;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
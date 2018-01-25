%%%%%%%%%%%%%%%%%%
% Rough confidence scores for reaction
% Reactions with pubmedID and with gene information -3
% Reactions with gene but without pubmedID -2; 
% Reactions without gene but need for modelling -1;
% Reactions without gene -0;
% Exchange reactions -none;
%%%%%%%%%%%%%%%%%%%%
function rxnConfidenceScores = getConfidenceScores(model)
t=length(model.rxnNames);
rxnConfidenceScores = cell(t,1);

for i = 1:t
    if ~isempty(model.rules{i}) 
        if ~isempty(strfind(model.rxnNotes{i},'pmid:'))
            rxnConfidenceScores{i} = 3;
        else
            rxnConfidenceScores{i} = 2;
        end
    else
            rxnConfidenceScores{i} = 0;
    
    end
        
end

for i = 1:t
    if ~isempty(strfind(model.rxnNames{i},'exchange'))
           rxnConfidenceScores{i} = nan;
    elseif ~isempty(strfind(model.rxnNames{i},'isa'))
           rxnConfidenceScores{i} = 1;
    else 
           rxnConfidenceScores{i} = rxnConfidenceScores{i};
    end
        
end






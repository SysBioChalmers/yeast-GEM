%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = addSBOterms(model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = addSBOterms(model)

%Get RAVEN model for matching names & compartments
model_r = ravenCobraWrapper(model);

%Add SBO terms for mets:
model.metSBOTerms = cell(size(model.mets));
for i = 1:length(model.mets)
    metName = model_r.metNames{i};
    if ismember(metName,{'biomass','DNA','RNA','protein','carbohydrate','lipid','cofactor','ion'}) ...
            || endsWith(metName,' backbone') || endsWith(metName,' chain')
        model.metSBOTerms{i} = 'SBO:0000649';     %Biomass
    else
        model.metSBOTerms{i} = 'SBO:0000247';     %Simple chemical
    end
end

%Add SBO terms for rxns:
model.rxnSBOTerms = cell(size(model.rxns));
for i = 1:length(model.rxns)
    rxnName   = model_r.rxnNames{i};
    metNames  = model_r.metNames(model.S(:,i) ~= 0);
    metComps  = model_r.metComps(model.S(:,i) ~= 0);
    metStoich = model_r.S(model.S(:,i) ~= 0,i);
    
    if length(metNames) == 1
        if strcmp(model_r.comps{metComps},'e')
            model.rxnSBOTerms{i} = 'SBO:0000627';	%Exchange rxn
            
        elseif metStoich > 0
            model.rxnSBOTerms{i} = 'SBO:0000628';	%Demand rxn
        else
            model.rxnSBOTerms{i} = 'SBO:0000632';	%Sink rxn
        end
        
    elseif strcmp(rxnName,'biomass pseudoreaction')
        model.rxnSBOTerms{i} = 'SBO:0000629';       %Biomass pseudo-rxn
        
    elseif strcmp(rxnName,'non-growth associated maintenance reaction')
        model.rxnSBOTerms{i} = 'SBO:0000630';       %ATP maintenance
        
    elseif contains(rxnName,'pseudoreaction') || contains(rxnName,'SLIME rxn')
        model.rxnSBOTerms{i} = 'SBO:0000395';       %Encapsulating process
        
    elseif length(unique(metComps)) > 1 && length(unique(metNames)) < length(metNames)
        model.rxnSBOTerms{i} = 'SBO:0000655';       %Transport rxn
        
    else
        model.rxnSBOTerms{i} = 'SBO:0000176';       %Metabolic rxn
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

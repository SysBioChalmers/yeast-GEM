% model = addSBOterms(model)
function model = addSBOterms(model)

%Define SBO terms for mets
metsSBO=cell(size(model.mets));
for i = 1:length(model.mets)
    metName = model.metNames{i};
    if ismember(metName,{'biomass','DNA','RNA','protein','carbohydrate','lipid','cofactor','ion'}) ...
            || endsWith(metName,' backbone') || endsWith(metName,' chain')
        metsSBO{i} = 'SBO:0000649';     %Biomass
    else
        metsSBO{i} = 'SBO:0000247';     %Simple chemical
    end
end

%Define SBO terms for rxns
rxnSBO = cell(size(model.rxns));
rxnSBO(:) = {'SBO:0000176'};       %Metabolic rxn, if nothing else
% Exchange, sink & demand (only 1 reactant)
reactantNumber=sum(model.S~=0,1);
reactantNumber=find(reactantNumber==1);
for i=1:numel(reactantNumber)
    idx=reactantNumber(i);
    if strcmp(model.comps{model.metComps(find(model.S(:,idx)))},'e') || ...
            strcmp(model.compNames{model.metComps(find(model.S(:,idx)))},'extracellular')
        rxnSBO{idx} = 'SBO:0000627'; %Exchange rxn
    elseif sum(model.S(:,idx))<0
        rxnSBO{idx} = 'SBO:0000632';	%Sink rxn
    else
        rxnSBO{idx} = 'SBO:0000628';	%Demand rxn
    end
end
% Transport reactions
i=getTransportRxns(model);
rxnSBO(i) = {'SBO:0000655'};
% Pseudo reactions
for i=numel(model.rxns)
    if strcmp(model.rxnNames(i),'biomass pseudoreaction')
        rxnSBO{i} = 'SBO:0000629';       %Biomass pseudo-rxn
    elseif strcmp(model.rxnNames(i),'non-growth associated maintenance reaction')
        rxnSBO{i} = 'SBO:0000630';       %ATP maintenance
    elseif contains(model.rxnNames(i),'pseudoreaction') || contains(model.rxnNames(i),'SLIME rxn')
        rxnSBO{i} = 'SBO:0000395';       %Encapsulating process
    end
end

%Define SBO terms for genes
geneSBO = cell(size(model.genes));
geneSBO(:) = {'SBO:0000243'};       %Gene

% Add SBO term if it wasn't annotated yet
model=editMiriam(model,'met','all','sbo',metsSBO,'fill');
model=editMiriam(model,'rxn','all','sbo',rxnSBO,'fill');
model=editMiriam(model,'gene','all','sbo',geneSBO,'fill');

end

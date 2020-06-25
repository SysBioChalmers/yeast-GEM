%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% checkMetBalance(model,metName,flux,show_all)
% Shows for a given metabolite all involved reactions on all compartments.
%
% INPUT:    model       A GEM as a structure object
%           metName     The name of the desired metabolite (e.g. 'NADH')
%           flux        (OPTIONAL) A solution of the model (nx1 vector)
%           show_all    (OPTIONAL) If false, only non-zero fluxes will be
%                       displayed (default = false; if no flux is given as
%                       input then default = true)
%
% OUTPUT:   On the command window displays all involved reactions,
%           including the flux value, the reaction ID and the formula. 
%           Results are displayed by compartment, and separated depending
%           if the metabolite is being consumed, produced or transported 
%           in each rxn. If no flux was given as input, a zero will appear
%           instead.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function checkMetBalance(model,metName,flux,show_all)

lb = model.lb;
ub = model.ub;

if nargin < 4
    show_all = false;
    if nargin < 3
        flux     = zeros(size(model.rxns));
        show_all = true;
    end
end

pos = zeros(size(model.mets));
for i = 1:length(model.metNames)
    if strcmp(model.metNames{i},metName)
        pos(i) = 1;
    end
end
pos = find(pos);

for i = 1:length(pos)
    rxns_cons = {};
    rxns_prod = {};
    for j = 1:length(model.rxns)
        if model.S(pos(i),j) ~= 0 && (show_all || abs(flux(j)) > 1e-10)
            rxn         = model.rxns{j};
            rxn_formula = printRxnFormula(model,model.rxns(j),false,false,true);
            rxn_formula = rxn_formula{1};
            if model.S(pos(i),j) < 0
                rxns_cons = [rxns_cons;{lb(j) flux(j) ub(j) rxn rxn_formula}];
            else
                rxns_prod = [rxns_prod;{lb(j) flux(j) ub(j) rxn rxn_formula}];
            end
        end
    end
    if ~isempty(rxns_cons) || ~isempty(rxns_prod)
        print_group(rxns_cons,'Consumption')
        print_group(rxns_prod,'Production')
    end
end

fprintf('\n')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_group(rxns,name)

if ~isempty(rxns)
    fprintf([name ':\n'])
    [~,order] = sort(cell2mat(rxns(:,1)),'descend');
    rxns      = rxns(order,:);
    for i = 1:length(order)
        fprintf([sprintf('%6.2e',rxns{i,1}) '\t' sprintf('%6.2e',rxns{i,2}) ...
            '\t' sprintf('%6.2e',rxns{i,3}) '\t' rxns{i,4} '\t\t' rxns{i,5} '\n'])
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
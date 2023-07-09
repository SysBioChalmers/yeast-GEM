%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = changeRules(model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = changeRules(model)

%Change and's & or's:
model.rules = strrep(model.rules,'&','and');
model.rules = strrep(model.rules,'|','or');

%Change gene ids:
for i = 1:length(model.genes)
    model.rules = strrep(model.rules,['x(' num2str(i) ')'],model.genes{i});
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
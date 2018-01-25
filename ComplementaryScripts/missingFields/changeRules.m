%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = changeRules(model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = changeRules(model)

for i = 1:length(model.rules)
    rule = model.rules{i};
    
    %Change and's & or's:
    rule = strrep(rule,'&','AND');
    rule = strrep(rule,'|','OR');
    
    %Change gene ids:
    start_pos = strfind(rule,'x(');
    while ~isempty(start_pos)
        end_pos   = strfind(rule,')');
        end_pos   = end_pos(end_pos > start_pos(1));
        gene_pos  = str2double(rule(start_pos(1)+2:end_pos(1)-1));
        if start_pos(1) == 0
            rule = [model.genes{gene_pos} rule(end_pos(1)+1:end)];
        elseif end_pos(1) == length(rule)
            rule = [rule(1:start_pos(1)-1) model.genes{gene_pos}];
        else
            rule = [rule(1:start_pos(1)-1) model.genes{gene_pos} rule(end_pos(1)+1:end)];
        end
        start_pos = strfind(rule,'x(');
    end
    
    model.rules{i} = rule;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
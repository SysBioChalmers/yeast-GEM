%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subSystem = findSubSystem(uniprot,kegg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subSystem = findSubSystem(uniprot,kegg)

subSystem = '';

%Gather all matches:
uniprot = strsplit(uniprot,' ');
for i = 1:length(uniprot)
    for j = 1:length(kegg{1})
        if strcmp(uniprot{i},kegg{1}{j}) && ~isempty(kegg{6}{j})
            subSystem = strjoin({subSystem,kegg{6}{j}},' ');
        end
    end
end

if ~isempty(subSystem)
    subSystem = subSystem(2:end);
    
    %Remove repetitions and weird stuff:
    subSystem    = strsplit(subSystem,' sce0');
    first_sub    = subSystem{1};
    subSystem{1} = first_sub(5:end);
    sub_array    = unique(subSystem);
    subSystem    = '';
    for i = 1:length(sub_array)
        sub_i = sub_array{i};
        weird = strcmp(sub_i,'0591  Linoleic acid metabolism') || ...
                strcmp(sub_i,'0590  Arachidonic acid metabolism') || ...
                strcmp(sub_i,'0592  alpha-Linolenic acid metabolism') || ...
                strcmp(sub_i,'0565  Ether lipid metabolism') || ...
                strcmp(sub_i,'0460  Cyanoamino acid metabolism') || ...
                strcmp(sub_i,'0680  Methane metabolism');
        if ~weird
            subSystem = strjoin({subSystem,sub_i},' ;; sce0');
        end
    end
    subSystem = subSystem(5:end);
end

if ~isempty(subSystem)
    subSystem = strsplit(subSystem,' ;; ');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [uni,EC] = findInDB(rxn_pos,model,DB)
% Matches the uniprot and EC number for a given rxn into a given database.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [uni,EC] = findInDB(rxn_pos,model,DB)

%Find simple gene sets for reaction:
[gene_sets,~] = getAllPath(model,model.rxns(rxn_pos));
uni           = cell(size(gene_sets));
EC            = cell(size(gene_sets));

for i = 1:length(gene_sets)
    %Split the gene set and match each gene:
    gene_sets{i} = strrep(gene_sets{i},'and','AND');
    gene_set = strsplit(gene_sets{i},' AND ');
    uni_set  = cell(size(gene_set));
    EC_set   = cell(size(gene_set));
    for j = 1:length(gene_set)
        for k = 1:length(DB{3})
            if sum(strcmp(DB{3}{k},gene_set{j})) > 0
                uni_set{j} = DB{1}{k};
                if ~isempty(DB{4}{k})
                    new_EC_set = strsplit(DB{4}{k},' ');
                    for l = 1:length(new_EC_set)
                        EC_set{j} = [EC_set{j} new_EC_set{l} ';'];
                    end
                end
            end
        end
        if isempty(EC_set{j})
            EC_set{j} = '';
        else
            EC_set{j} = EC_set{j}(1:end-1);
        end
    end
    %Uniprot: Delete repeated and empty spaces
    [uni_set,~] = deleteRepeated(uni_set);
    uni_set     = uni_set(~cellfun('isempty',uni_set));
    
    %EC: Find union and intersection between all units (only applies for
    %complexes, i.e. length(EC_set) > 1):
    uni_EC = strsplit(EC_set{1},' ');
    int_EC = uni_EC;
    for j = 2:length(EC_set)
        other_EC = strsplit(EC_set{j},' ');
        if isempty(uni_EC)
            uni_EC = other_EC;
            int_EC = other_EC;
        elseif ~isempty(other_EC)
            uni_EC = compare_wild([uni_EC other_EC]);
            int_EC = intersection(int_EC,other_EC);
        end
    end
    %Use the intersection found, if any. If not, use the union: 
    if isempty(int_EC)
        EC_set = uni_EC;
    else
        EC_set = int_EC;
    end
        
    %Add new codes as new possible isoenzymes:
    uni{i} = union_string(uni_set);
    EC{i}  = union_string(EC_set);
end

%Create array with all options:
uni = union_string(uni);
EC  = union_string(EC);
uni = strsplit(uni,' ');
EC  = strsplit(EC,';');

%Delete repeated stuff:
uni = unique(uni);
EC  = unique(EC);

%Return everything as a string:
uni = union_string(uni);
EC  = union_string(EC);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function int_EC = intersection(prev_EC,new_EC)
%Finds the common elements between two cell arrays, if any. Also considers
%wildcards (e.g. if 'EC1.1.1.1' is in one array and 'EC1.1.1.-' is in the
%other one, then 'EC1.1.1.1' is added to the intersection).

int_EC = {};
for i = 1:length(prev_EC)
    for j = 1:length(new_EC)
        new_int = compare_wild({prev_EC{i} new_EC{j}});
        if length(new_int) == 1
            int_EC = [int_EC new_int];
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function int_EC = compare_wild(EC)
%Goes through a cell array of EC numbers, and erases any repetitions,
%considering also wildcards (e.g. will erase 'EC1.2.3.-' if 'EC1.2.3.4' is
%already present).

%Trim all EC numbers of wild cards (e.g. 'EC1.2.-.-' -> 'EC1.2.'):
EC_trimmed = EC;
for i = 1:length(EC)
    %Add a final dot for avoiding issues (e.g. 'EC1.1.1.1' & 'EC1.1.1.12'):
    ECi = [EC{i} '.'];
    pos = strfind(ECi,'-');
    if ~isempty(pos)
        ECi = ECi(1:pos(1)-1);
    end
    EC_trimmed{i} = ECi;
end

%Compare all EC numbers between them to find if 1 fits in the other:
non_repeated = true(1,length(EC));
for i = 1:length(EC)-1
    for j = i+1:length(EC)
        ECi = EC_trimmed{i};
        ECj = EC_trimmed{j};
        %If ECj fits in ECi then ECj can be disregarded:
        if contains(ECi,ECj)
            non_repeated(j) = false;
        %Else, if ECi fits in ECj then ECi can be disregarded:
        elseif contains(ECj,ECi)
            non_repeated(i) = false;
        end
    end
end

int_EC = EC(non_repeated);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function str = union_string(cell_array)
%Receives any 1xn cell array and returns the union of all non empty
%elements as a string

nonempty = ~cellfun(@isempty,cell_array);
str      = strjoin(cell_array(nonempty)',';');

end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
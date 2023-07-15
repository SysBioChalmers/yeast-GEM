%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [list,deleted] = deleteRepeated(list)
% Recieves a list, removes all repeated elements and returns the list with
% only non repeated elements, and also the list of the positions deleted.
%
% NOTE: Must be a string or cell array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [list,deleted] = deleteRepeated(list)

N       = length(list);
deleted = false(1,N);
for i = 1:N-1
    for j = i+1:N
        if strcmp(list{i},list{j})
            deleted(j) = true;
        end
    end
end
list(deleted) = [];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
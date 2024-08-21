function model = nitrogenLimitation(model)
% nitrogenLimitation
%   Converts model to represents nitrogen-limiting environmental conditions
%
%   Inputs: model           (struct) unmodified model
%   Output: model           (struct) nitrogen-limitation model
%
%   Usage: model = nitrogenLimitation(model)

% Glutamine synthase is repressed when nitrogen is in excess. See doi:10.1128/EC.2.5.827-829.2003
model.ub(strcmp(model.rxns,'r_0472'))=1000;
end

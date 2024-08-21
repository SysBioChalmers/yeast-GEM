function model = glycineNitrogenSource(model)
% glycineNitrogenSource
%   Converts model to represent glycine as sole nitrogen source
%
%   Inputs: model           (struct) unmodified model
%   Output: model           (struct) glycine model
%
%   Usage: model = glycineNitrogenSource(model)

% Glycine cleavage is only active when glycine is used as sole nitrogen source. See doi:10.1111/j.1567-1364.2002.tb00069.x and doi:10.1074/jbc.274.15.10523 
doi:10.1128/EC.2.5.827-829.2003
model.ub(strcmp(model.rxns,'r_0501'))=0; %glycine cleavage, mitochondrion
model.lb(strcmp(model.rxns,'r_0501'))=0;
model.ub(strcmp(model.rxns,'r_0507'))=0; %glycine cleavage complex (lipoylprotein), mitochondrion
model.lb(strcmp(model.rxns,'r_0507'))=0;
model.ub(strcmp(model.rxns,'r_0509'))=0; %glycine cleavage complex (lipoamide), mitochondrion
model.lb(strcmp(model.rxns,'r_0509'))=0;
end
    
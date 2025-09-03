function model = scaleBioMass(model,component,new_value,balance_out,dispOutput)
% scaleBioMass
%   Scales the biomass composition
%
% Input:
%   model           (struct) yeast-GEM model
%   component       (string) biomass component to change (options are:
%                   'carbohydrate', 'protein', 'lipid', 'RNA', 'DNA',
%                   'ion', 'cofactor')
%   new_value       (num) new total fraction for the specified biomass
%                   component
%   balance_out     (string, optional) biomass component that will be used
%                   to balance out the biomass composition, so that the
%                   total mass adds up to 1 g/gDCW. This is highly
%                   recommended (default = empty, no scaling takes place)
%   dispOutput      (bool, optional) displayed outoupt (default = true)
%
% Output:
%   model          (struct) modified yeast-GEM model
%
% Usage: model = scaleBioMass(model,component,new_value,balance_out,dispOutput)

if nargin < 5
    dispOutput = true;
end
if nargin < 4
    balance_out = '';
end
  
%Measure current composition and rescale:
[~,P,C,R,D,L,I,F] = sumBioMass(model,dispOutput);
content_all = {'carbohydrate','protein','lipid','RNA','DNA','ion','cofactor'};
content_Cap = {'C','P','L','R','D','I','F'};
pos         = strcmp(content_all,component);
old_value   = eval(content_Cap{pos});
f           = new_value / old_value;
model       = rescalePseudoReaction(model,component,f);
X           = sumBioMass(model,false);

%Balance out (if desired):
if ~isempty(balance_out)
    pos           = strcmp(content_all,balance_out);
    balance_value = eval(content_Cap{pos});
    f             = (balance_value + (1-X)) / balance_value;
    model         = rescalePseudoReaction(model,balance_out,f);
end
end

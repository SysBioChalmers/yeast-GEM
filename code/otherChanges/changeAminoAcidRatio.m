function model = changeAminoAcidRatio(model,aerobic)
% changeAminoAcidRatio
%   Updates the amino acid ratio in the biomass equation, based on data by
%   Björkeroth et al. (2020, PNAS, doi:10.1073/pnas.1921890117). 
%
% Input:
%   model   (struct)  the yeast GEM
%   aerobic (logical) true if ratios should be set for aerobic condition,
%                     false for anaerobic conditions. Default is true.
%
% Output:
%   model   (struct)  the updated yeast GEM
%
% Usage: model = changeAminoAcidRatio(model,aerobic)
%

if nargin < 2 || isempty(aerobic)
    col = 1;
elseif aerobic % First column is aerobic, second column anaerobic
    col = 1;
else
    col = 2;
end

funcDir = dbstack('-completenames');
funcDir = regexprep(funcDir(1).file,[funcDir(1).name '\.m'],'');

%Load chemostat data:
fid     = fopen(fullfile(funcDir,'..','..','data','physiology','aminoacid_Bjorkeroth2020.tsv'),'r');
data    = textscan(fid,'%s %s %s %f %f %f','Delimiter','\t','HeaderLines',1);
tRNAids = [data{2} data{3}];
aaRatio = [data{5} data{6}];
fclose(fid);

aaRatio = aaRatio(:,col);
aaRatio = [-aaRatio; aaRatio];

[~, P] = sumBioMass(model,false);
protRxn  = find(strcmp(model.rxns,'r_4047'));
[~, tRNAidxs] = ismember(tRNAids,model.mets);
model.S(tRNAidxs,protRxn) = aaRatio;
model = scaleBioMass(model,'protein',P,[],false);
end
function model = changeAminoAcidRatio(model,col)

if nargin < 2 || isempty(col)
    col = 1;
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

protRxn  = getIndexes(model,'r_4047','rxns');
[~, tRNAidxs] = ismember(tRNAids,model.mets);
model.S(tRNAidxs,protRxn) = aaRatio;

model = scaleBioMass(model,'protein',P,[],false);
end
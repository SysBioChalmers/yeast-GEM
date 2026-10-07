function dumpEnvironments(outDir)
% dumpEnvironments  Applies every environment in data/conditions with
%   applyEnvironment and writes, per environment, <name>_bounds.tsv (rxn, lb,
%   ub) and <name>_S.tsv (met, rxn, coef of every nonzero). The Python side
%   compares these with yeastgem.conditions.apply:
%   python compare_environments.py <outDir>
%
%   Needs yeast-GEM's code folder and a toolbox that loads model/yeast-GEM.yml
%   (loadYeastYaml) on the path; applyEnvironment itself needs neither.

if ~isfolder(outDir), mkdir(outDir); end
codeDir = fileparts(which('applyEnvironment'));
model = loadYeastYaml;
files = dir(fullfile(codeDir,'..','data','conditions','*.yml'));
for i = 1:numel(files)
    name = erase(files(i).name,'.yml');
    m = applyEnvironment(model,name);
    writetable(table(m.rxns,m.lb,m.ub,'VariableNames',{'rxn','lb','ub'}), ...
        fullfile(outDir,[name '_bounds.tsv']),'FileType','text','Delimiter','\t');
    [r,c,v] = find(m.S);
    writetable(table(m.mets(r),m.rxns(c),full(v),'VariableNames',{'met','rxn','coef'}), ...
        fullfile(outDir,[name '_S.tsv']),'FileType','text','Delimiter','\t');
    fprintf('%s written\n',name);
end
end

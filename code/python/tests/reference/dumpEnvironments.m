function dumpEnvironments(outDir)
% dumpEnvironments  Applies every environment in data/conditions with
%   applyEnvironment and writes, per environment, <name>.tsv in outDir
%   (default: the folder environments next to this file): the reactions whose
%   bounds changed (kind bounds, id, lb, ub) and the stoichiometric
%   coefficients that changed (kind coef, metabolite, reaction, new value).
%   test_environments_match_matlab.py applies the same environments with
%   yeastgem.conditions.apply and requires the same changes.
%
%   The changes are relative to the model with every exchange reaction opened
%   (lb -1000, ub 1000).
%
%   Needs yeast-GEM's code folder and a toolbox that loads model/yeast-GEM.yml
%   (loadYeastYaml) on the path; applyEnvironment itself needs neither.
%   Run it again whenever a file in data/conditions changes.

if nargin < 1
    outDir = fullfile(fileparts(mfilename('fullpath')),'environments');
end
if ~isfolder(outDir), mkdir(outDir); end
codeDir = fileparts(which('applyEnvironment'));
model = loadYeastYaml;
% Open every exchange reaction first: the model file already holds the minimal
% medium, so environments that reset the exchanges would change nothing.
isExchange = ~any(model.S > 0,1)' | ~any(model.S < 0,1)';
model.lb(isExchange) = -1000;
model.ub(isExchange) = 1000;
files = dir(fullfile(codeDir,'..','data','conditions','*.yml'));
for i = 1:numel(files)
    name = erase(files(i).name,'.yml');
    m = applyEnvironment(model,name);
    lines = {};
    for j = find(m.lb ~= model.lb | m.ub ~= model.ub)'
        lines{end+1,1} = sprintf('bounds\t%s\t\t%.12g\t%.12g',m.rxns{j},m.lb(j),m.ub(j)); %#ok<AGROW>
    end
    [r,c] = find(abs(m.S - model.S) > 1e-12);
    for j = 1:numel(r)
        lines{end+1,1} = sprintf('coef\t%s\t%s\t%.12g\t',m.mets{r(j)},m.rxns{c(j)}, ...
            full(m.S(r(j),c(j)))); %#ok<AGROW>
    end
    fid = fopen(fullfile(outDir,[name '.tsv']),'w');
    fprintf(fid,'kind\tid\tid2\tvalue\tvalue2\n');
    lines = sort(lines);
    fprintf(fid,'%s\n',lines{:});
    fclose(fid);
    fprintf('%s: %d changes\n',name,numel(lines));
end
end

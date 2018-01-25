%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saveYeastModel(model)
% Saves model as a .mat, .xml and .txt file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function saveYeastModel(model)

model.description = 'yeastGEM';

%Remove any space in rxnECNumbers:
model.rxnECNumbers = strrep(model.rxnECNumbers,' ','');

%Save changes to current model:
writeCbModel(model,'sbml','../ModelFiles/xml/yeastGEM.xml');
writeCbModel(model,'text','../ModelFiles/txt/yeastGEM.txt');

%Detect boundary metabolites and save them in a .txt file:
fid = fopen('boundaryMets.txt','wt');
for i = 1:length(model.rxns)
    pos = find(model.S(:,i) ~= 0);
    if length(pos) == 1 %Exchange rxn
        fprintf(fid,[model.mets{pos} '\t' model.metNames{pos} '\n']);
    end
end
fclose(fid);

%Retrieve SBML toolbox version:
SBMLpath = which('SBMLToolbox.m');
slashPos = strfind(SBMLpath,'\');
if isempty(slashPos)
    slashPos = strfind(SBMLpath,'/');
end
try
    SBMLpath = SBMLpath(1:slashPos(end-1));
    fid      = fopen([SBMLpath 'VERSION.txt'],'r');
    SBMLTver = fscanf(fid,'%s');
    fclose(fid);
catch
    SBMLTver = '?';
end

%Save file with versions:
fid = fopen('dependencies.txt','wt');
fprintf(fid,['SBML_toolbox\tv' SBMLTver '\n']);
fields = fieldnames(model.modelVersion);
for i = 1:length(fields)
    value = model.modelVersion.(fields{i});
    fprintf(fid,[fields{i} '\t' num2str(value) '\n']);
end
fclose(fid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
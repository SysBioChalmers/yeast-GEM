%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saveYeastModel(model,version)
% Saves model as an .sbml and .txt file.
%
% Benjamín J. Sánchez. Last edited: 2017-10-31
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function saveYeastModel(model,version)

model.description = ['yeast_' version '.xml'];

%Remove any space in rxnECNumbers:
model.rxnECNumbers = strrep(model.rxnECNumbers,' ','');

%Save changes to current model:
cd ../ModelFiles/mat
save(['yeast_' version '.mat'],'model');
cd ../xml
writeCbModel(model,'sbml',['yeast_' version '.xml']);
cd ../txt
writeCbModel(model,'text',['yeast_' version '.txt']);

%Detect boundary metabolites and save them in a .txt file:
cd ../../ComplementaryScripts
fid = fopen('boundaryMets.txt','wt');
for i = 1:length(model.rxns)
    pos = find(model.S(:,i) ~= 0);
    if length(pos) == 1 %Exchange rxn
        fprintf(fid,[model.mets{pos} '\t' model.metNames{pos} '\n']);
    end
end
fclose(fid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
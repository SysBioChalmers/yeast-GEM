%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% increaseVersion(model,version)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function increaseVersion(model,version)

%Include tag in model:
model.description = ['yeastGEM_v' version];
saveYeastModel(model)

%Store model as .mat (only for releases):
mkdir('../ModelFiles/mat');
save('../ModelFiles/mat/yeastGEM.mat','model');

%Update version file:
fid = fopen('../version.txt','wt');
fprintf(fid,version);
fclose(fid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
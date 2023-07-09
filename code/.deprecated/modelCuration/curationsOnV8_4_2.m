%% Curations on version 8.4.2
% This is a list of various smaller curations to be performed on version 
% 8.4.2 of yeast-GEM. Indicated is what Issue is solved.

% Load model
cd ..
model = loadYeastModel;

% Solves #252
model = changeGeneAssociation(model,'r_4590','YOL130W');

% Solves #254
model.metNames = regexprep(model.metNames,'^nicotinamide ribose','nicotinamide riboside');
% All annotations, formula and charge were correct.

saveYeastModel(model);

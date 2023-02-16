cd ..
seed = readcell("seeddata.xlsx");
dgpredictor = readcell("rxn_G.xlsx");
cd model_from_yetfl
load("yeast8_thermo_curated.mat");
cd ../../code
yeast = loadYeastModel;
tmpmetNames = yeast.metNames;
for i = 1:length(yeast.mets)
    tmpmetNames{i, 1} = strcat(yeast.metNames{i, 1}, ' [', yeast.compNames{yeast.metComps(i, 1), 1}, ']');
end
yeast.metsdeltaG = {};
% if the mets in yetfl, choose the G in yetfl (https://github.com/EPFL-LCSB/yetfl)
for i = 1:length(model.mets)
    yeast.metsdeltaG{i, 1} = model.metDeltaGFtr(find(strcmp(tmpmetNames{i, 1}, model.metNames)));
end
% if the mets not in yetfl, choose the G in ModelSeed database
for i = length(model.mets) + 1:length(yeast.mets)
    kegg = find(strcmp('kegg.compound', yeast.metMiriams{i, 1}.name));
    if kegg
        g = seed{find(strcmp(yeast.metMiriams{i, 1}.value{kegg,1}, seed(:, 6))), 4};
        yeast.metsdeltaG{i, 1} = str2double(g);
    end
end

% if the rxn in yetfl, choose the G in yetfl (https://github.com/EPFL-LCSB/yetfl)
for i = 1:3982
    yeast.rxnsdeltaG{i, 1} = model.rxnDeltaGR(find(strcmp(yeast.rxns{i, 1}, model.rxns)));
end

% if the rxn not in yetfl, choose the G calculated by dgpredictor (https://doi.org/10.1371/journal.pcbi.1009448)
for i = 3982:length(yeast.rxns)
    if length(dgpredictor{i+1, 3}) == 1
        yeast.rxnsdeltaG{i, 1} = dgpredictor{i+1, 3};
    end
end
cd ../deltaG/yeast-GEM_with_G
save yeast;


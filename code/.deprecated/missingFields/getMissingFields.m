%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = getMissingFields(model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = getMissingFields(model)

%Change rule format:
old_rules = model.rules;
model = changeRules(model);

%Load swissprot:
fid       = fopen('../../data/databases/swissprot.tsv','r');
swissprot = textscan(fid,'%s %s %s %s %f32 %s','Delimiter','\t','HeaderLines',1);
swissprot = standardizeDatabase(swissprot);
fclose(fid);

%Load kegg:
fid  = fopen('../../data/databases/kegg.tsv','r');
kegg = textscan(fid,'%s %s %s %s %f32 %s %s','Delimiter','\t','HeaderLines',1);
kegg = standardizeDatabase(kegg);
fclose(fid);

%Get rxn EC numbers and uniprots:
model.rxnECNumbers = cell(size(model.rxns));
model.subSystems   = cell(size(model.rxns));
classification     = zeros(length(model.rxns),7);
for i = 1:length(model.rxns)
    %Find matches for UNIPROT codes and EC numbers:
    [new_uni_swiss,new_EC_swiss] = findInDB(i,model,swissprot);
    [new_uni_kegg,new_EC_kegg]   = findInDB(i,model,kegg);
    
    %Asign EC numbers (prioritize Swissprot):
    if isempty(new_EC_swiss)
        model.rxnECNumbers{i} = new_EC_kegg;
    else
        model.rxnECNumbers{i} = new_EC_swiss;
    end
    
    %Assign subSystem (prioritize KEGG):
    if ~isempty(new_uni_kegg)
        model.subSystems{i} = findSubSystem(new_uni_kegg,kegg);
    end
    
    if isempty(model.subSystems{i}) && ~isempty(new_uni_swiss)
        model.subSystems{i} = findSubSystem(new_uni_swiss,kegg);
    end
    
    %Clasify reaction:
    %1st column: is it an exchange rxn?
    classification(i,1) = sum(model.S(:,i)~=0)==1;
    %2nd column: is it a transport rxn?
    met_pos = model.mets(model.S(:,i)~=0);
    for j = 1:length(met_pos)-1
        met_j  = met_pos{j};
        pos_j  = strfind(met_j,'[');
        comp_j = met_j(pos_j(end)+1:end-1);
        for k = j+1:length(met_pos)
            met_k  = met_pos{k};
            pos_k  = strfind(met_k,'[');
            comp_k = met_k(pos_k(end)+1:end-1);
            if ~strcmp(comp_j,comp_k)
                classification(i,2) = 1;
            end
        end
    end
    %3rd column: does the rxn have gene associations?
    classification(i,3) = ~isempty(model.rules{i});
    %4th column: does the rxn have EC number?
    classification(i,4) = ~isempty(model.rxnECNumbers{i});
    %5th column: is the rxn part of a subsystem?
    classification(i,5) = ~isempty(model.subSystems{i});
    %6th column: more than 1 EC number?
    classification(i,6) = length(strfind(model.rxnECNumbers{i},';'))+1;
    %7th column: more than 1 subsystem?
    classification(i,7) = length(strfind(model.subSystems{i},'sce0'));
    disp(['Adding missing fields: Ready with rxn ' int2str(i)])
end

%Show main counts:
tot = length(model.rxns);
disp(['Total of rxns: ' num2str(tot)])
exc = sum(classification(:,1));
disp(['Total of exchange rxns: ' num2str(exc) ' (' num2str(exc/tot*100) '% of total)'])
tra = sum(classification(:,2));
disp(['Total of transport rxns: ' num2str(tra) ' (' num2str(tra/tot*100) '% of total)'])
enz = sum(sum(classification(:,1:2),2)==0);
disp(['Total of enzymatic rxns: ' num2str(enz) ' (' num2str(enz/tot*100) '% of total)'])
gen = sum(classification(:,3));
disp(['Total of rxns with gene info: ' num2str(gen) ' (' num2str(gen/tot*100) '% of total)'])
ecn = sum(classification(:,4));
disp(['Total of rxns with EC number: ' num2str(ecn) ' (' num2str(ecn/gen*100) '% coverage)'])
sub = sum(classification(:,5));
disp(['Total of rxns with subsystem: ' num2str(sub) ' (' num2str(sub/gen*100) '% coverage)'])
mec = sum(classification(:,6)>1);
disp(['Rxns with >1 EC numbers: ' num2str(mec) ' (' num2str(mec/ecn*100) '% of all)'])
msu = sum(classification(:,7)>1);
disp(['Rxns with >1 subsystems: ' num2str(msu) ' (' num2str(msu/sub*100) '% of all)'])

%Save model (with old rules):
cd ..
model.rules = old_rules;
saveYeastModel(model)
cd missingFields

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function database = standardizeDatabase(database)

for i = 1:length(database{3})
    database{3}{i} = strsplit(database{3}{i},' ');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
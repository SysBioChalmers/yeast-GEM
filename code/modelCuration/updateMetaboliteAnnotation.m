%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = updateMetaboliteAnnotation(model)
% update the metabolite annotation information in the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = updateMetaboliteAnnotation(model)

%Load data:
fid = fopen('../../data/modelCuration/metabolite_manual_curation.tsv','r');
metaboliteData = textscan(fid,'%s %s %s %s %f32 %s','Delimiter','\t','HeaderLines',1);
fclose(fid);

for i = 1:length(metaboliteData{1})
    for j = 1:length(model.mets)
        %Correct name:
        if startsWith(model.metNames{j},[metaboliteData{1}{i} ' ['])	%old name
            metName = model.metNames{j};
            comp    = metName(strfind(metName,' ['):end);
            metName = [metaboliteData{2}{i} comp];
            model.metNames{j} = metName;
        end
        
        %Update other fields:
        if startsWith(model.metNames{j},[metaboliteData{2}{i} ' ['])	%new name
            model.metChEBIID{j} = metaboliteData{3}{i};	%new CHEBI
            model.metKEGGID{j}  = metaboliteData{4}{i};	%new KEGG
            model.metCharges(j) = metaboliteData{5}(i);	%new charge
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
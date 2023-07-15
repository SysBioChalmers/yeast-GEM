%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convertYmn2FBC2
% The main purpose of this script is to convert yeast metabolic network 7.6
% from COBRA-compatible SBML2 to FBC v2 thereby adding the missing
% annotation data, which could not be retained with the older COBRA
% versions. The annotation covers metabolites and reactions.
% So we cannot save gene annotation data for now. There are also problems
% with KEGG reaction ids and reaction literature references, see below.
% We also had to remove "-" characters from gene names to "_", since
% otherwise gene-reaction association information is lost for several
% reactions in FBC v2 file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Firstly loading an old COBRA version, i.e. the latest commit before
% pulling request from 'aebrahim-patch-1' branch, which brings FBC v2
% support and thereby discontinues SBML2 support.
% So the last SBML2-compatible COBRA version can be downloaded from
% https://github.com/opencobra/cobratoolbox/tree/10f98a0ab834118c6dfb943970c38d68c7e1ae70
% 
% The whole idea of the script is to load curated Yeast 7.6 from 'main'
% branch, then load SBML3 format model from
% https://sourceforge.net/projects/yeast/files/, add the missing
% annotation to the former model, do some final changes and export it to
% FBC v2.
model_sbml2=readCbModel('yeast_7.6_cobra.xml');
save('yeast_7.6_cobra_sbml2.mat','model_sbml2');

% Now resetting Matlab pathlist and adding up-to-date COBRA to path
matlabrc

% Re-opening the model again
load('yeast_7.6_cobra_sbml2.mat');
% Now we open Yeast 7.6 model in SBML3 (FBC v1) format, obtained from sourceforge
model_sbml3=readCbModel('yeast_7.6.xml');

% Now let's check whether reactions and metabolites are in the same order
isequal(model_sbml2.rxns,model_sbml3.rxns)
sbml3_mets = regexprep(model_sbml3.mets,'\[.+\]',''); % removing [..] tailing;
isequal(model_sbml2.mets,sbml3_mets)
isequal(model_sbml2.grRules,model_sbml3.grRules)
isequal(model_sbml2.rxnGeneMat,model_sbml3.rxnGeneMat)
% Everything is fine, so now we can append SBML2 model with the data
% obtained from SBML3.
clear ans sbml3_mets

% We can already copy metabolite annotation data, which was imported by
% readCbModel function
model_sbml2.metCHEBIID=model_sbml3.metCHEBIID;
model_sbml2.metKEGGID=model_sbml3.metKEGGID;
model_sbml2.metPubChemID=model_sbml3.metPubChemID;
model_sbml2.metInChIString=model_sbml3.metInChIString;

% We no longer need model_sbml3;
clear model_sbml3

% The following information was not imported by readCbModel:
%   model.metCharge. We will import it.
%   KEGG reaction ids. There is no designated variable in model structure
%   for it, so we cannot carry it.
%   model.rxnReferences. When having multiple references for one reaction,
%   separated by semicolon in FBC v2, only the first reference is imported,
%   when using readCbModel. So it is better for now to leave it empty
%   rather than having it incomplete.

% As for metCharge, we need to open SBML3 model with TranslateSBML function
% and to manually obtain this data. modelSBML.species contains metabolites
% and genes, so we will need to map modelSBML.species to model_sbml2.mets
modelSBML=TranslateSBML('yeast_7.6.xml');
metCharge=cell(numel(model_sbml2.mets),1);
for i=1:numel(model_sbml2.mets)
    for j=1:numel(modelSBML.species)
        if strcmp(model_sbml2.mets(i),modelSBML.species(j).id)
            metCharge{i}=regexprep(modelSBML.species(j).notes,'.+<p>CHARGE:','');
            metCharge{i}=regexprep(metCharge{i},'</p>.+','');
            metCharge{i}=regexprep(metCharge{i},'+','');
        end
    end
end
model_sbml2.metCharge=str2double(metCharge);
clear i j metCharge

% Now let's obtain KEGG reaction ids and literature references for
% reactions
% rxnReferences=cell(numel(model_sbml2.rxns),1);

% UPDATE: The following loop was written to fetch KEGG reaction ids and literature
% references. It is not possible for now to keep this annotation data, but
% the loop may still be useful later, if aforementioned issues are fixed.
%
%for i=1:numel(modelSBML.reaction)
%	rxnReferences{i}=regexprep(modelSBML.reaction(i).annotation,'.+http://identifiers.org/kegg.reaction/','http://identifiers.org/kegg.reaction/');
%    rxnReferences{i}=regexprep(rxnReferences{i},'"/>         </rdf:Bag>       </bqbiol:is>       <bqbiol:isDescribedBy>         <rdf:Bag>           <rdf:li rdf:resource="','; ');
%    rxnReferences{i}=regexprep(rxnReferences{i},'"/>           <rdf:li rdf:resource="','; ');
%    % Now let's assume that particular reaction doesn't have KEGG id;
%    rxnReferences{i}=regexprep(rxnReferences{i},'.+ <rdf:li rdf:resource="http://identifiers.org/pubmed/','http://identifiers.org/pubmed/');
%    rxnReferences{i}=regexprep(rxnReferences{i},'"/>.+','');
%end
%
%model_sbml2.rxnReferences=rxnReferences;
%clear i modelSBML rxnReferences;


% We are almost ready to export model to FBC v2 format

% For successful exporting model.mets must be in
% e.g. H2O[c] format. We obtain compartments from model.metNames
compartments=regexprep(model_sbml2.metNames,'(.+)\[','');
compartments=strcat('[',compartments,'');

model_sbml2.mets=strcat(model_sbml2.mets,compartments,'');
model_sbml2.mets=regexprep(model_sbml2.mets,'cell envelope','ce'); % cell envelope
model_sbml2.mets=regexprep(model_sbml2.mets,'cytoplasm','c'); % cytoplasm
model_sbml2.mets=regexprep(model_sbml2.mets,'endoplasmic reticulum membrane','erm'); % e.r. membrane
model_sbml2.mets=regexprep(model_sbml2.mets,'endoplasmic reticulum','er'); % e.r.
model_sbml2.mets=regexprep(model_sbml2.mets,'extracellular','e'); % extracellular
model_sbml2.mets=regexprep(model_sbml2.mets,'Golgi membrane','gm'); % Golgi membrane
model_sbml2.mets=regexprep(model_sbml2.mets,'Golgi','g'); % Golgi
model_sbml2.mets=regexprep(model_sbml2.mets,'lipid particle','lp'); % lipid particle
model_sbml2.mets=regexprep(model_sbml2.mets,'mitochondrial membrane','mm'); % mitochondrial membrane
model_sbml2.mets=regexprep(model_sbml2.mets,'mitochondrion','m'); % mitochondria
model_sbml2.mets=regexprep(model_sbml2.mets,'nucleus','n'); % nucleus
model_sbml2.mets=regexprep(model_sbml2.mets,'peroxisome','p'); % peroxisome
model_sbml2.mets=regexprep(model_sbml2.mets,'vacuolar membrane','vm'); % vacuolar membrane
model_sbml2.mets=regexprep(model_sbml2.mets,'vacuole','v'); % vacuole

model_sbml2=rmfield(model_sbml2,'modelVersion');

% Now fixing gene names;
model_sbml2.genes=regexprep(model_sbml2.genes,'-','_');
model_sbml2.grRules=regexprep(model_sbml2.grRules,'-','_');

% Save model:
saveYeastModel(model_sbml2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
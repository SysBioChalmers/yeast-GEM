% This script is adds new database annotated reactions to the model.
% =MORE DETAILED DESCRIPTION REQUIRED. WHICH DATABASES? WHICH DATE? WHAT
% CRITERIA OR SEARCH TERMS WERE USED? ANY ADDITIONAL FILTERING?=
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd ..
model = loadYeastModel;
metsInfo = '../data/modelCuration/DBnewRxns/DBnewRxnsMets.tsv';
rxnsCoeffs = '../data/modelCuration/DBnewRxns/DBnewRxnsCoeffs.tsv';
rxnsInfo = '../data/modelCuration/DBnewRxns/DBnewRxnsRxns.tsv';
genesInfo = '../data/modelCuration/DBnewRxns/DBnewRxnsGenes.tsv';
newModel = curateMetsRxnsGenes(model,metsInfo,genesInfo,rxnsCoeffs,rxnsInfo);
checkModelStruct(newModel,true,false)
saveYeastModel(newModel)
cd modelCuration
newModel=model;
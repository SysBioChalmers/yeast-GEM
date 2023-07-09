% This script adds missing metabolites/reactions involving sulfur metabolism,
% as described in Huang et al. 2017 FEMS YR doi:10.1093/femsyr/fox058
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd ..
model = loadYeastModel;
metsInfo = '../data/modelCuration/volatileSulfurCompounds/VoSulMets.tsv';
rxnsCoeffs = '../data/modelCuration/volatileSulfurCompounds/VoSulRxnsCoeffs.tsv';
rxnsInfo = '../data/modelCuration/volatileSulfurCompounds/VoSulRxns.tsv';
genesInfo = '../data/modelCuration/volatileSulfurCompounds/VoSulGenes.tsv';
newModel = curateMetsRxnsGenes(model,metsInfo,genesInfo,rxnsCoeffs,rxnsInfo);
saveYeastModel(model)
cd modelCuration

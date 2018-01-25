%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Ptot,Ctot] = calculateContent(model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ptot,Ctot] = calculateContent(model)

% MW aminoacids [g/mol]:
aas = {'s_0404[c]'	89.09       % A     Alanine         ala
       's_0542[c]'	121.16      % C     Cysteine        cys
       's_0432[c]'	133.11      % D     Aspartic acid   asp
       's_0748[c]'	147.13      % E     Glutamic acid   glu
       's_1314[c]'	165.19      % F     Phenylalanine   phe
       's_0757[c]'	75.07       % G     Glycine         gly
       's_0832[c]'	155.15      % H     Histidine       his
       's_0847[c]'	131.17      % I     Isoleucine      ile
       's_1099[c]'	146.19      % K     Lysine          lys
       's_1077[c]'	131.17      % L     Leucine         leu
       's_1148[c]'	149.21      % M     Methionine      met
       's_0430[c]'	132.12      % N     Asparagine      asn
       's_1379[c]'	115.13      % P     Proline         pro
       's_0747[c]'	146.14      % Q     Glutamine       gln
       's_0428[c]'	174.2       % R     Arginine        arg
       's_1428[c]'	105.09      % S     Serine          ser
       's_1491[c]'	119.12      % T     Threonine       thr
       's_1561[c]'	117.15      % V     Valine          val
       's_1527[c]'	204.23      % W     Tryptophan      trp
       's_1533[c]'	181.19};    % Y     Tyrosine        tyr

% MW carbohidrates [g/mol]:
carbs = {'s_0001[ce]'	180.16      % (1->3)-beta-D-glucan
         's_0004[ce]'	180.16      % (1->6)-beta-D-glucan
         's_0509[c]'	221.21      % chitin
         's_0773[c]'	180.16      % glycogen
         's_1107[c]'	180.16      % mannan
         's_1520[c]'	342.296};	% trehalose

%Initialize protein and carb content:
Ptot = 0;
Ctot = 0;

%Count protein/carb content in the corresponding pseudo-rxn:
protPos = strcmp(model.metNames,'protein [cytoplasm]');
carbPos = strcmp(model.metNames,'carbohydrate [cytoplasm]');
protRxn = model.S(protPos,:) == 1;
carbRxn = model.S(carbPos,:) == 1;
for i = 1:length(model.mets)
    posP = strcmp(aas(:,1),model.mets{i});
    posC = strcmp(carbs(:,1),model.mets{i});
    if sum(posP) == 1
        Sprot = abs(model.S(i,protRxn));            % mmol/gDW
        Ptot  = Ptot + Sprot*aas{posP,2}/1000;      % mmol/gDW * g/mmol = g/gDW
    elseif sum(posC) == 1
        Scarb = abs(model.S(i,carbRxn));            % mmol/gDW
        Ctot  = Ctot + Scarb*carbs{posC,2}/1000;    % mmol/gDW * g/mmol = g/gDW
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
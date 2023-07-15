function rxn_lnRevIdx = get_lnRevIdx(deltaG,rxn,model)
% get_lnRevIdx
%   Calculates the log-scale reversibility index of a reaction
%   Reference to paper by Noor et al., doi:10.1093/bioinformatics/bts317
%   
%   deltaG      Gibbs free energy of reaction
%   rxn         Identifier in model.rxns e.g. r_0001
%   model       COBRA model structure
%
%   Usage: rxn_lnRevIdx = get_lnRevIdx(deltaG,rxn,model)
%

%define constants
R = 8.3145e-3; %gas constant (in kJ/K·mol)
T = 298.15; %temperature at standard conditions (in K)
C = 1e-3; %characteristic physiological concentration (in M)
if nargin~=3
    warning('Incorrect input: deltaG and corresponding identifier model.rxns required');
else
    if isa(deltaG,'double')
        deltaG = deltaG*4.184; %unit conversion from kcal/mol to kJ/mol
        %find number of reactants
        rxn_idx = find(ismember(model.rxns,rxn));
        met_idx = find(model.S(:,rxn_idx));
        coef = model.S(met_idx,rxn_idx);
        N_sub = 0;
        N_pdt = 0;
        for j = 1:length(coef)
            if coef(j) <= -1
                N_sub = N_sub + 1;
            elseif coef(j) >= 1
                N_pdt = N_pdt + 1;
            end
        end
        %calculate reversibility index (in natural logarithm scale)
        N = N_pdt + N_sub;
        if (N_pdt - N_sub) ~= 0
            %calculate physiological deltaG correction
            physio_deltaG_correction = (N_pdt - N_sub)*log(C);
            rxn_lnRevIdx = ((2/N)*(deltaG) + R*T*physio_deltaG_correction)/(R*T);
        elseif (N_pdt - N_sub) == 0
            rxn_lnRevIdx = ((2/N)*(deltaG))/(R*T);
        end
    end
end

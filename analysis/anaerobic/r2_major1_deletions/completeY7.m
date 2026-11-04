function model = completeY7(model)
% completeY7  Copy of the complete_Y7 local function of yeast-GEM v9.1.0
%   code/modelTests/essentialGenes.m: Kennedy synthetic complete medium.
    % change Y7 model medium to the Kennedy synthetic complete medium
    % NOTE: Y7 doesn't have exchange reactions for dihydrofolate,
    % octadecenoate, octadecynoate, or linoleic acid.

    % start with a clean slate: set all exchange reactions to upper bound =
    % 1000 and lower bound = 0 (ie, unconstrained excretion, no uptake)
    [~, exchangeRxns] = getExchangeRxns(model);
    model.lb(exchangeRxns) = 0;
    model.ub(exchangeRxns) = 1000;

    constrainedUptake = {'r_1604';'r_1639';'r_1873';'r_1879';'r_1880';...
        'r_1881';'r_1671';'r_1883';'r_1757';'r_1891';'r_1889';'r_1810';...
        'r_1993';'r_1893';'r_1897';'r_1947';'r_1899';'r_1900';'r_1902';...
        'r_1967';...
        %'r_1918';... Y7 doesn't have a linoleic acid exchange rxn, which
        %was substituted for octadecenoate and octadecynoate for Y5 and Y6,
        %so leave out.
        'r_1903';'r_1548';'r_1904';'r_2028';...
        'r_2038';'r_1906';'r_2067';'r_1911';'r_1912';'r_1913';'r_2090';...
        'r_1914';'r_2106'};   % potassium exchange};
    glucoseExchange = {'r_1714';};
    unconstrainedUptake = {'r_1672';'r_1654'; ... % ammonium exchange
                    'r_1992'; ... % oxygen exchange
                    'r_2005'; ... % phosphate exchange
                    'r_2060'; ... % sulphate exchange
                    'r_1861'; ... % iron exchange, for test of expanded biomass def
                    'r_1832'; ... % hydrogen exchange
                    'r_2100'; ... % water exchange
                    'r_4593'; ... % chloride exchange
                    'r_4595'; ... % Mn(2+) exchange
                    'r_4596'; ... % Zn(2+) exchange
                    'r_4597'; ... % Mg(2+) exchange
                    'r_2049'; ... % sodium exchange
                    'r_4594'; ... % Cu(2+) exchange
                    'r_4600'; ... % Ca(2+) exchange
                    'r_2020' };

    model=setParam(model,'lb',constrainedUptake,-0.5);
    model=setParam(model,'lb',glucoseExchange,-20);
    model=setParam(model,'lb',unconstrainedUptake,-1000);
end


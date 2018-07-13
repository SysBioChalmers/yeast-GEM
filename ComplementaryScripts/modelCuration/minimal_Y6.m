function model = minimal_Y6(model)
    % change Y6 model media to minimal - ammonium, glucose, oxygen,
    % phosphate, sulphate
    % the function is from:https://doi.org/10.1371/journal.pcbi.1004530
    % start with a clean slate: set all exchange reactions to
    % upper bound = 1000 and lower bound = 0 (ie, unconstrained
    % excretion, no uptake)
    exchangeRxns = findExcRxns(model);
    model.lb(exchangeRxns) = 0;
    model.ub(exchangeRxns) = 1000;
    desiredExchanges = {...
        'r_1654'; ... % 'ammonium exchange';
        'r_1992'; ... % 'oxygen exchange';
        'r_2005'; ... % 'phosphate exchange';
        'r_2060'; ... % 'sulphate exchange';
        'r_1861'; ... % iron for test of expanded biomass def;
        'r_1832'; ... % hydrogen exchange;
        };
    glucoseExchange = {...
        'r_1714'; ... % D-glucose exchange'
        };
    uptakeRxnIndexes = findRxnIDs(model,desiredExchanges);
    glucoseExchangeIndex = findRxnIDs(model,glucoseExchange);
    if length(uptakeRxnIndexes) ~= 6 %5;
        error('Not all exchange reactions were found.')
    end
    model.lb(uptakeRxnIndexes)=-1000;
    model.lb(glucoseExchangeIndex)=-10;
end
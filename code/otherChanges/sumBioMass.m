function [X,P,C,R,D,L,I,F] = sumBioMass(model,dispOutput)
  % sumBioMass
  %   Calculates breakdown of biomass
  %
  %   model         (struct) Metabolic model in COBRA format
  %   dispOutput    (bool, opt) If output should be displayed (default = true)
  %
  %   X             (float) Total biomass fraction [gDW/gDW]
  %   P             (float) Protein fraction [g/gDW]
  %   C             (float) Carbohydrate fraction [g/gDW]
  %   R             (float) RNA fraction [g/gDW]
  %   D             (float) DNA fraction [g/gDW]
  %   L             (float) Lipid fraction [g/gDW]
  %   F             (float) cofactor [g/gDW]
  %   I             (float) ion [g/gDW]
  %
  %   Usage: [X,P,C,R,D,L,I,F] = sumBioMass(model,dispOutput)
  %
  %   Function adapted from SLIMEr: https://github.com/SysBioChalmers/SLIMEr
  %

if nargin < 2
    dispOutput = true;
end

%Get main fractions:
[P,X] = getFraction(model,'P',0,dispOutput);
[C,X] = getFraction(model,'C',X,dispOutput);
[R,X] = getFraction(model,'R',X,dispOutput);
[D,X] = getFraction(model,'D',X,dispOutput);
[L,X] = getFraction(model,'L',X,dispOutput);
[I,X] = getFraction(model,'I',X,dispOutput);
[F,X] = getFraction(model,'F',X,dispOutput);

if dispOutput
    disp(['X -> ' num2str(X) ' gDW/gDW'])
    % Simulate growth:
    sol = solveLP(model,1);
    disp(['Growth = ' num2str(sol.f) ' 1/h'])
    disp(' ')
end
end

%%
function [F,X] = getFraction(model,compType,X,dispOutput)

%Define pseudoreaction name:
rxnName = [compType ' pseudoreaction'];
rxnName = strrep(rxnName,'P','protein');
rxnName = strrep(rxnName,'C','carbohydrate');
rxnName = strrep(rxnName,'N','biomass');
rxnName = strrep(rxnName,'L','lipid backbone');
rxnName = strrep(rxnName,'R','RNA');
rxnName = strrep(rxnName,'D','DNA');
rxnName = strrep(rxnName,'I','ion');
rxnName = strrep(rxnName,'F','cofactor');

%Add up fraction:
rxnPos = strcmp(model.rxnNames,rxnName);
if isempty(rxnPos)
    if dispOutput
        disp([compType ' does not exist '])
    end
    F = 0;
else
    isSub   = find(model.S(:,rxnPos)<0); % Substrates in pseudoreaction
    if strcmp(compType,'L') % Lipid already has g/gDW as unit
        F = full(-sum(model.S(isSub,rxnPos)));
    else
        formulas = model.metFormulas(isSub);
        MWs = zeros(numel(formulas),1);
        for i = 1:numel(formulas)
            MWs(i) = parseChemicalFormula(formulas{i});
        end
        zeroMW = MWs == 0;
        if any(zeroMW)
            error('Biomass metabolite %s has an empty metFormula field.', model.mets{isSub(zeroMW)})
        end
        switch compType
            case 'P'
                % Two protons have to be removed from the charged-tRNA
                % formulas that are in the model
                MWs = MWs - 2.016;
            case {'R','D'}
                % H2O has to be removed to represent polymerization
                MWs = MWs - 18.015;
        end
        F = full(-sum(model.S(isSub,rxnPos).*MWs)/1000);
    end
end
X = X + F;

if dispOutput
    disp([compType ' -> ' num2str(F) ' g/gDW'])
end
end

function molecularWeight = parseChemicalFormula(formula)
    % Split formula in elements and coefficients
    tokens = regexp(formula, '([A-Z][a-z]*)(\d*)', 'tokens');
    tokensMatrix = vertcat(tokens{:});
    tokensMatrix(cellfun(@isempty,tokensMatrix(:,2)),2) = {'1'};
    elements = tokensMatrix(:, 1);
    counts   = str2double(tokensMatrix(:, 2));

    %Weight of elements
    elem    = {'C', 12.01; ...
               'H', 1.008; ...
               'N', 14.007; ...
               'O', 15.999; ...
               'P', 30.974; ...
               'S', 32.06; ...
               'R', 0; ...
               'Fe', 55.845; ...
               'K', 39.098; ...
               'Na', 22.99; ...
               'Cl', 35.45; ...
               'Mn', 54.938; ...
               'Zn', 65.38; ...
               'Ca', 40.078; ...
               'Mg', 24.305; ...
               'Cu', 63.546};
    [~,elemMatch] = ismember(elements,elem(:,1));
    molecularWeight = sum(counts .* transpose([elem{elemMatch,2}]),'all');
end

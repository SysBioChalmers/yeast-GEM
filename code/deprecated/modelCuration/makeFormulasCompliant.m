%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makeFormulasCompliant.m
% Corrects SBML-incompatible metabolite formulas.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%aa codes:
aas = {'s_0955[c]'	'ala'   	% A     Alanine
       's_0981[c]'	'cys'       % C     Cysteine
       's_0973[c]'	'asp'       % D     Aspartic acid
       's_0991[c]'	'glu'       % E     Glutamic acid
       's_1032[c]'	'phe'   	% F     Phenylalanine
       's_1003[c]'	'gly'       % G     Glycine
       's_1006[c]'	'his'   	% H     Histidine
       's_1016[c]'	'ile'       % I     Isoleucine
       's_1025[c]'	'lys'       % K     Lysine
       's_1021[c]'	'leu'   	% L     Leucine
       's_1029[c]'	'met'       % M     Methionine
       's_0969[c]'	'asn'   	% N     Asparagine
       's_1035[c]'	'pro'       % P     Proline
       's_0999[c]'	'gln'   	% Q     Glutamine
       's_0965[c]'	'arg'       % R     Arginine
       's_1039[c]'	'ser'   	% S     Serine
       's_1045[c]'	'thr'       % T     Threonine
       's_1056[c]'	'val'       % V     Valine
       's_1048[c]'	'trp'       % W     Tryptophan
       's_1051[c]'	'tyr'};     % Y     Tyrosine

%tRNA's are just produced and consumed in single loops, so we can just
%replace the radical part with an "R" and still maintain mass balances,
%i.e. aa-tRNA(aa) will have "R" + the formula of the corresponding aa and
%tRNA(aa) will just have an "R". Example:
%     tRNA(Ala): C10H17 O10PR2(C5H8O6PR)n -> R
% Ala-tRNA(Ala): C13H22NO11PR2(C5H8O6PR)n -> C3H6NOR (ala-R without an -OH)
% Cycle in which the 2 are involved:
% r_4047: 0.4193 Ala-tRNA(Ala) + ...  -> 0.4193 tRNA(Ala) + ... + protein (growth)
% r_0157: ATP           + L-alanine + tRNA(Ala) -> Ala-tRNA(Ala) + AMP         + diphosphate
%         C10H12N5O13P3   C3H7NO2     R      	   C3H6NOR		   C10H12N5O7P	 HO7P2

for i = 1:length(model.mets)
    name = model.metNames{i};
    if contains(name,'-tRNA(')
        %Correct metabolite:
        try
            aaName = lower(name(1:strfind(name,'-')-1));
            aaID   = aas{strcmp(aas(:,2),aaName),1};
            aaPos  = strcmp(model.mets,aaID);
            model.metFormulas{i} = [model.metFormulas{aaPos} 'R'];
            model.metFormulas{i} = takeOutFromFormula(model.metFormulas{i},'OH');
        catch
            %formyl-met: C6H11NO3S
            model.metFormulas{i} = 'C6H10NO2SR';
        end
        
        %Correct associated tRNA pair:
        pairName = name(strfind(name,'-')+1:end);
        pairPos  = strcmp(model.metNames,pairName);
        if sum(pairPos) > 0
            model.metFormulas{pairPos} = 'R';
        end
    end
end

%Try to use COBRA's "computeMetFormulae":
problems = contains(model.metFormulas,'(');
model_test = model;
model_test.metFormulas(problems) = '';
[~,metFormulae,~,~,~,~,~,~] = computeMetFormulae(model_test);

for i = 1:length(model.mets)
    if problems(i)
        name        = model.metNames{i};
        formula     = model.metFormulas{i};
        formula_new = '';
        rxn_set     = find(model.S(i,:) ~= 0);
        disp(' ')
        disp([model.mets{i} ' - ' name ':']);
        try
            pos = strcmp(metFormulae(:,1),model.mets{i});
            formula_new = metFormulae{pos,2};
        catch
            disp('Could not balance with COBRA!')
        end
        
        %For the rest solve by manual inspection of rxn formulas:
        for j = rxn_set
            %Substrates:
            rxn_subs = model.S(:,j) < 0;
            rxn_subs = model.metFormulas(rxn_subs);
            rxn_subs = strjoin(rxn_subs,' + ');
            %Products:
            rxn_prods = model.S(:,j) > 0;
            rxn_prods = model.metFormulas(rxn_prods);
            rxn_prods = strjoin(rxn_prods,' + ');
            disp([model.rxns{j} ' - ' model.rxnNames{j} ': ' rxn_subs ' -> ' rxn_prods])
        end
        
        %Solve depending on the case:
        if strcmp(formula_new,'')
            if contains(name,'-beta-D-glucan ')
                %All glucans can have the same formula:
                formula_new = 'C6H10O5';
            elseif contains(name,'dolichol ') || contains(name,'dolichyl ')
                %Remove "(C5H8)n", as it doesn't respect mass balances:
                formula_new = formula(1:strfind(formula,'(')-1);
            elseif contains(name,'chitin ') || contains(name,'chitosan ') || ...
                    contains(name,'pectin ')
                %Keep inner part, as the rest is "H2O(...)" to simbolyze that is 
                %hydrated, but doesn't respect mass balances:
                formula_new = formula(strfind(formula,'(')+1:strfind(formula,')')-1);
            elseif contains(name,'lipoylprotein ')
                %Remove "(C2H2NOR)n", as it doesn't respect mass balances:
                formula_new = formula(1:strfind(formula,'(')-1);
            end
        end
        
        %Adapt metabolite formula:
        model.metFormulas{i} = formula_new;
        disp(['Resolution: ' formula ' -> ' formula_new]);
        
        %Verify that reactions associated to the metabolite are balanced with RAVEN:
        for j = rxn_set
            getElementalBalance(model,j,true,true);
        end
        pause
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
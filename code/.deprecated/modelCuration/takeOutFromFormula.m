%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% formula = takeOutFromFormula(formula,elements)
% Takes away from formula each of the elements specified.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function formula = takeOutFromFormula(formula,elements)

for i = 1:length(elements)
    posInFormula = strfind(formula,elements(i));
    if posInFormula == length(formula)
        formula = strrep(formula,elements(i),'');               %only 1 element at the end of the formula
    else
        for j = posInFormula+1:length(formula)
            if isletter(formula(j))
                if j == posInFormula+1
                    formula = strrep(formula,elements(i),'');   %only 1 element in the middle of the formula
                else
                    formula = goDown(formula,posInFormula,j);	%more than 1 element in the formula
                end
                break
            elseif j == length(formula)
                formula = goDown(formula,posInFormula,j+1);     %more than 1 element at the end of the formula
            end
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function formula = goDown(formula,startPos,endPos)

numbElem = str2double(formula(startPos+1:endPos-1));
if numbElem == 2
    numbElem = '';
else
    numbElem = num2str(numbElem - 1);
end
if endPos == length(formula)+1
    formula = [formula(1:startPos) numbElem];
else
   formula = [formula(1:startPos) numbElem formula(endPos:end)];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
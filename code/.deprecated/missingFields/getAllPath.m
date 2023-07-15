%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [FinalGenesets, FinalReactions] = getAllPath(model,reaction)
% Generates clean gene sets with only simple relations (enzyme alone or
% complexes, no isoenzymes). Used by getECnumbers for distinguishing
% complex GPR relationships.
%
% INPUT:
% model             The GEM structure (1x1 struct)
% reaction          The reaction ID   (string)
%
% OUTPUTS:
% FinalGenesets     All gene sets with only "ANDs" (cell array)
% FinalReactions    New reaction IDs (cell array)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [FinalGenesets, FinalReactions] = getAllPath(model,reaction)
i = ismember(model.rxns,reaction);        %get index of the selected reaction

FinalGenesets = cell(1,1);
FinalReactions = cell(1,1);

x = 0;      %total number of NewGenesets
NewGenesets = cell(1,1);
NewReactions = cell(1,1);
STR = model.rules{i};
p = 0;
q = 0;
while p < x+1      %if bracket exist in geneset, means there still exist multi-gene controlled reaction, continue iteration
    STR = strtrim(STR);
    GRbracketL = strfind(STR,'(');     %find all indexes of left bracket
    GRbracketR = strfind(STR,')');     %find all indexes of right bracket
    A = zeros(length(GRbracketL),1);    %ready for saving corresponding dual index of left brackets
    B = zeros(length(GRbracketR),1);    %ready for saving corresponding dual index of right brackets    
         %save the geneset as string STR
    tempSTR = STR;
    for k = 1:length(A)     %identification of all dual brackets
        B(k,1) = min(GRbracketR);   %iteratively get the index of corresponding dual right brackets
        A(k,1) = max(strfind(tempSTR(1:min(GRbracketR)),'('));  %iteratively get the index of corresponding dual left brackets
        GRbracketR = setdiff(GRbracketR,min(GRbracketR));
        tempSTR = [tempSTR(1:A(k,1)-1) '[' tempSTR(A(k,1)+1:length(tempSTR))];
    end
    if  min(A) == 1 & B(A(:,1)==min(A)) == length(STR)     %if there are overall brackets
        STR = STR(min(A)+1:B(A(:,1)==min(A))-1);     %remove the overall bracket
    end
    STR = strtrim(STR);
    GRbracketLtemp = strfind(STR,'(');    %remove the index of overal left bracket
    GRbracketRtemp = strfind(STR,')');    %remove the index of overal right bracket
    A = zeros(length(GRbracketLtemp),1);    %ready for saving corresponding dual index of left brackets
    B = zeros(length(GRbracketRtemp),1);    %ready for saving corresponding dual index of right brackets
    TempGenesets = cell(1,1);
    if ~isempty(A)      %if there are inner brackets 
            
        tempSTR = STR;
        for k = 1:length(A)     %identification of all dual brackets
            B(k,1) = min(GRbracketRtemp);   %iteratively get the index of corresponding dual right brackets
            A(k,1) = max(strfind(tempSTR(1:min(GRbracketRtemp)),'('));  %iteratively get the index of corresponding dual left brackets
            GRbracketRtemp = setdiff(GRbracketRtemp,min(GRbracketRtemp));
            tempSTR = [tempSTR(1:A(k,1)-1) '[' tempSTR(A(k,1)+1:length(tempSTR))];
        end
        subSTR = STR;        
        j = 0;       
        while ~isempty(strfind(subSTR,'('))                
                j = j+1;
                TempGenesets{j,1} = subSTR(min(A):B(A(:,1)==min(A))); %add new Geneset to Geneset
               
                if min(A) <= 1
                    if B(A(:,1)==min(A)) < length(subSTR)
                        subSTR = subSTR(B(A(:,1)==min(A))+1:length(subSTR));
                    end
                elseif B(A(:,1)==min(A)) >= length(subSTR)
                    if min(A) > 1
                        subSTR = subSTR(1:min(A)-1);
                    end
                else
                    subSTR = [subSTR(1:min(A)-1) subSTR(B(A(:,1)==min(A))+1:length(subSTR))];
                end
                
                GRbracketLtemp2 = strfind(subSTR,'(');    %remove the index of overal left bracket
                GRbracketRtemp2 = strfind(subSTR,')');    %remove the index of overal right bracket
                A = zeros(length(GRbracketLtemp2),1);    
                B = zeros(length(GRbracketRtemp2),1);
                tempSTR2 = subSTR;
                for k = 1:length(A)     %identification of all dual brackets in subSTR
                    B(k,1) = min(GRbracketRtemp2);   %iteratively get the index of corresponding dual right brackets
                    A(k,1) = max(strfind(tempSTR2(1:min(GRbracketRtemp2)),'('));  %iteratively get the index of corresponding dual left brackets
                    GRbracketRtemp2 = setdiff(GRbracketRtemp2,min(GRbracketRtemp2));
                    tempSTR2 = [tempSTR2(1:A(k,1)-1) '[' tempSTR2(A(k,1)+1:length(tempSTR2))];
                end
        end             
        for m = 1:length(model.genes)
            if ~isempty(strfind(subSTR,model.genes{m}))               
                TempGenesets{length(TempGenesets)+1,1} = model.genes{m};
            end    
        end 
        if ~isempty(strfind(subSTR,' or ')) || ~isempty(strfind(subSTR,' OR '))
            for k = 1:length(TempGenesets)
                x = x+1;
                NewGenesets{x,1} = TempGenesets{k,1};
                NewReactions{x,1} = {['TempReactions' num2str(x)]};
            end
        elseif ~isempty(strfind(subSTR,' and ')) || ~isempty(strfind(subSTR,' AND '))
            n = 0;
            jumpout = 0;
            genesetSTR = STR;
            while n < j && jumpout == 0
                Left = strfind(genesetSTR,'(');
                Right = strfind(genesetSTR,')');
                A3 = zeros(length(Left),1);    %ready for saving corresponding dual index of left brackets
                B3 = zeros(length(Right),1);
                tempgenesetSTR = genesetSTR;
                for k = 1:length(A3)     %identification of all dual brackets
                    B3(k,1) = min(Right);   %iteratively get the index of corresponding dual right brackets
                    A3(k,1) = max(strfind(tempgenesetSTR(1:min(Right)),'('));  %iteratively get the index of corresponding dual left brackets
                    Right = setdiff(Right,min(Right));
                    tempgenesetSTR = [tempgenesetSTR(1:A3(k,1)-1) '[' tempgenesetSTR(A3(k,1)+1:length(tempgenesetSTR))];
                end        
                LeftBrackets = A3;
                RightBrackets = B3;                
                n = n+1;
                STR2 = TempGenesets{n};
                STR2 = strtrim(STR2);
                GRbracketL2 = strfind(STR2,'(');     %find all indexes of left bracket
                GRbracketR2 = strfind(STR2,')');     %find all indexes of right bracket               
                A2 = zeros(length(GRbracketL2),1);    %ready for saving corresponding dual index of left brackets
                B2 = zeros(length(GRbracketR2),1);    %ready for saving corresponding dual index of right brackets    
                tempSTR2 = STR2;
                for k = 1:length(A2)     %identification of all dual brackets
                    B2(k,1) = min(GRbracketR2);   %iteratively get the index of corresponding dual right brackets
                    A2(k,1) = max(strfind(tempSTR2(1:min(GRbracketR2)),'('));  %iteratively get the index of corresponding dual left brackets
                    GRbracketR2 = setdiff(GRbracketR2,min(GRbracketR2));
                    tempSTR2 = [tempSTR2(1:A2(k,1)-1) '[' tempSTR2(A2(k,1)+1:length(tempSTR2))];
                end
                if  min(A2) == 1 & B2(A2(:,1)==min(A2)) == length(STR2)
                    STR2 = STR2(min(A2)+1:B2(A2(:,1)==min(A2))-1);     %remove the overall bracket
                end  
                GRbracketLtemp2 = strfind(STR2,'(');    %remove the index of overal left bracket
                GRbracketRtemp2 = strfind(STR2,')');    %remove the index of overal right bracket
                A2 = zeros(length(GRbracketLtemp2),1);    %ready for saving corresponding dual index of left brackets
                B2 = zeros(length(GRbracketRtemp2),1);    %ready for saving corresponding dual index of right brackets
                TempGenesets2 = cell(1,1);
                if ~isempty(A2)      %if there are inner brackets             
                    tempSTR2 = STR2;
                    for k = 1:length(A2)     %identification of all dual brackets
                        B2(k,1) = min(GRbracketRtemp2);   %iteratively get the index of corresponding dual right brackets
                        A2(k,1) = max(strfind(tempSTR2(1:min(GRbracketRtemp2)),'('));  %iteratively get the index of corresponding dual left brackets
                        GRbracketRtemp2 = setdiff(GRbracketRtemp2,min(GRbracketRtemp2));
                        tempSTR2 = [tempSTR2(1:A2(k,1)-1) '[' tempSTR2(A2(k,1)+1:length(tempSTR2))];
                    end
                    subSTR2 = STR2;        
                    j2 = 0;      %recorder of gene sets
                    while ~isempty(strfind(subSTR2,'('))                
                        j2 = j2+1;
                        TempGenesets2{j2,1} = subSTR2(min(A2):B2(A2(:,1)==min(A2))); %add new Geneset to Geneset                
                        if min(A2) <= 1
                            if B2(A2(:,1)==min(A2)) < length(subSTR2)
                                subSTR2 = subSTR2(B2(A2(:,1)==min(A2))+1:length(subSTR2));
                            end
                        elseif B2(A2(:,1)==min(A2)) >= length(subSTR2)
                            if min(A2) > 1
                                subSTR2 = subSTR2(1:min(A2)-1);
                            end
                        else
                            subSTR2 = [subSTR2(1:min(A2)-1) subSTR2(B2(A2(:,1)==min(A2))+1:length(subSTR2))];
                        end
                        GRbracketLtemp3 = strfind(subSTR2,'(');    %remove the index of overal left bracket
                        GRbracketRtemp3 = strfind(subSTR2,')');    %remove the index of overal right bracket
                        A2 = zeros(length(GRbracketLtemp3),1);    
                        B2 = zeros(length(GRbracketRtemp3),1);
                        tempSTR2 = subSTR2;
                        for k = 1:length(A2)     %identification of all dual brackets in subSTR
                            B2(k,1) = min(GRbracketRtemp3);   %iteratively get the index of corresponding dual right brackets
                            A2(k,1) = max(strfind(tempSTR2(1:min(GRbracketRtemp3)),'('));  %iteratively get the index of corresponding dual left brackets
                            GRbracketRtemp3 = setdiff(GRbracketRtemp3,min(GRbracketRtemp3));
                            tempSTR2 = [tempSTR2(1:A2(k,1)-1) '[' tempSTR2(A2(k,1)+1:length(tempSTR2))];
                        end
                    end
                    for m = 1:length(model.genes)
                        if ~isempty(strfind(subSTR2,model.genes{m}))                
                            TempGenesets2{length(TempGenesets2)+1,1} = model.genes{m};                
                        end    
                    end
                else
                    subSTR2 = STR2;
                    j2 = 0;
                    for m = 1:length(model.genes)
                        if ~isempty(strfind(subSTR2,model.genes{m}))
                            TempGenesets2{j2+1,1} = model.genes{m};
                            j2 = j2+1;                                            
                        end    
                    end
                end
                if ~isempty(strfind(subSTR2,' or '))|| ~isempty(strfind(subSTR2,' OR '))
                    jumpout = 1;
                    for k = 1:length(TempGenesets2)
                        x = x+1;
                        if min(LeftBrackets) <= 1
                            if RightBrackets(LeftBrackets(:,1)==min(LeftBrackets)) < length(genesetSTR)
                                NewGenesets{x,1} = [TempGenesets2{k,1} genesetSTR(RightBrackets(LeftBrackets(:,1)==min(LeftBrackets))+1:length(genesetSTR))];
                            end
                        elseif RightBrackets(LeftBrackets(:,1)==min(LeftBrackets)) >= length(genesetSTR)
                            if min(LeftBrackets) > 1
                                NewGenesets{x,1} = [genesetSTR(1:min(LeftBrackets)-1) TempGenesets2{k,1}];
                            end
                        else
                            NewGenesets{x,1} = [genesetSTR(1:min(LeftBrackets)-1) TempGenesets2{k,1} genesetSTR(RightBrackets(LeftBrackets(:,1)==min(LeftBrackets))+1:length(genesetSTR))];
                        end
                        NewReactions{x,1} = {['TempReactions' num2str(x)]};                            
                    end
                else
                    if min(LeftBrackets) <= 1
                        if RightBrackets(LeftBrackets(:,1)==min(LeftBrackets)) < length(genesetSTR)
                            genesetSTR = [STR2 genesetSTR(RightBrackets(LeftBrackets(:,1)==min(LeftBrackets))+1:length(genesetSTR))];
                        end
                    elseif RightBrackets(LeftBrackets(:,1)==min(LeftBrackets)) >= length(genesetSTR)
                        if min(LeftBrackets) > 1
                            genesetSTR = [genesetSTR(1:min(LeftBrackets)-1) STR2];
                        end
                    else
                        genesetSTR = [genesetSTR(1:min(LeftBrackets)-1) STR2 genesetSTR(RightBrackets(LeftBrackets(:,1)==min(LeftBrackets))+1:length(genesetSTR))];
                    end    
                end
            end
            if jumpout == 0
                x = x+1;
                NewGenesets{x,1} = genesetSTR;
            end            
        end
        A = zeros(length(GRbracketLtemp),1); 
    elseif isempty(A)
        subSTR = STR;        
        if ~isempty(strfind(subSTR,' or '))|| ~isempty(strfind(subSTR,' OR '))
            n = 0;
            for m = 1:length(model.genes)
                if ~isempty(strfind(subSTR,model.genes{m}))
                    n = n+1;                
                    TempGenesets{n,1} = model.genes{m};                
                end    
            end 
            for k = 1:length(TempGenesets)
                q = q+1;
                FinalGenesets{q,1} = TempGenesets{k,1};
                FinalReactions{q,1} = {[reaction 'No' num2str(q)]};
            end
        elseif ~isempty(strfind(subSTR,' and '))|| ~isempty(strfind(subSTR,' AND '))
            q = q+1;
            FinalGenesets{q,1} = STR;
            FinalReactions{q,1} = {[reaction 'No' num2str(q)]};
        else
            q = q+1;
            FinalGenesets{q,1} = STR;
            FinalReactions{q,1} = {[reaction 'No' num2str(q)]};
        end
    end
    p = p+1;
    if p <= length(NewGenesets)
        STR = NewGenesets{p,1};
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
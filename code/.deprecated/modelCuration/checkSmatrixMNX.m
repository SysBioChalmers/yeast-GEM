% This script does an automated check of reaction balance in the model by
% comparing the S matrix of the model reactions
% to the S matrix of the MNX reactions via the MNXRID and MNXMID
%
% If a reaction is only missing either H+ or H2O (i.e. MNXM1 or MNXM2), and
% if the reaction in MNX database is balanced, then model.S will be
% modified accordingly to balance the equation
%
% Inputs: model, reac_prop.tsv
% Outputs: noMNXRID, noCoef, coefNoMatch, missingrxnmets, balanced and
% changeMNXM
%
% The outputs of this script can be used to facilitate manual curation of
% reaction balance in the model
%
% Note: This automatic curation had little impact in reactifying 
% reaction balance in the model, but the outputs were helpful 
% in the manual curation process of reaction balance
% The script was thus added into repository for documentation purposes
%

%Load model
cd ..
model = loadYeastModel;

%Check for reac_prop.tsv in current directory
%If not available, file will be downloaded
downloadMNXdb('reac_prop',pwd)

%Load reac_prop.tsv file
fileDir = dir('reac_prop.tsv');
if ~isempty(fileDir) %check if reactions.tsv is present in current directory
    fid = fopen('reac_prop.tsv');
    format = repmat('%s ',1,6);
    format = strtrim(format);
    rxn_temp = textscan(fid,format,'Delimiter','\t','HeaderLines',0);
    for i = 1:length(rxn_temp)
        MNXreac_prop(:,i) = rxn_temp{i};
    end
    commentLines = startsWith(MNXreac_prop(:,1),'#');
    MNXreac_prop(commentLines,:) = []; %MNXRID Equation Description Balance EC Source
    fclose(fid);
    
    %% Generate S matrix for MNX and model for comparison
    
    %Preallocate cell arrays
    noMNXRID{3991,2} = [];
    noCoef{3991,2} = [];
    coefNoMatch{3991,2} = [];
    allmetsDiff{3991,2} = [];
    balanced{3991,2} = [];
    changeMNXM{3991,2} = [];
    missingrxnmets{3991,2} = [];
    extrarxnmets{3991,2} = [];
    uncertainComp{3991,2} = [];
    modelR = ravenCobraWrapper(model);
    
    for i = 1:length(model.rxns)
        [~,rxnIdx] = ismember(model.rxnMetaNetXID(i),MNXreac_prop(:,1));
        if isequal(rxnIdx,0)
            arrayidx = find(cellfun('isempty', noMNXRID),1);
            noMNXRID(arrayidx,:) = [model.rxns(i), 'MNXRID not found, please check and modify to a new ID'];
        else
            rxn = strrep(MNXreac_prop(rxnIdx,2),'=','<=>');
            [S, mets]=constructS(rxn);
            matrix.metcoef = S;
            if isequal(matrix.metcoef,0)
                arrayidx = find(cellfun('isempty', noCoef),1);
                noCoef(arrayidx,:) = [model.rxns(i), 'Coefficient of mets are not found, please check and modify accordingly'];
            else
                temp = split(mets,'@',2);
                matrix.metMetaNetXID = cellstr(temp(:,1));
                
                %Fix the issue for H+: two IDs to represent H+ in MNX database
                matrix.metMetaNetXID(ismember(matrix.metMetaNetXID,'MNXM01')) = {'MNXM1'}; %
                MNXrefrxn = [matrix.metMetaNetXID,num2cell(matrix.metcoef)];
                coef = model.S(:,i);
                mets_temp = (find(model.S(:,i)));
                mets_MNX  = model.metMetaNetXID(mets_temp);
                mets_coef = coef(mets_temp);
                metComps = modelR.metComps;
                metComps = modelR.comps(metComps);
                mets_comp = metComps(mets_temp);
                modelrefrxn = [mets_MNX,num2cell(mets_coef),mets_comp];
                
                if any(ismember(mets_MNX,''))
                    arrayidx = find(cellfun('isempty', missingrxnmets),1);
                    missingrxnmets(arrayidx,:) = [model.rxns(i), 'mets_MNX contains empty cell array'];
                    
                elseif~all(ismember(mets_MNX,matrix.metMetaNetXID)) || ~all(ismember(matrix.metMetaNetXID,mets_MNX))
                    diffmets = setdiff(matrix.metMetaNetXID,mets_MNX);
                    diffmets2 = setdiff(mets_MNX, matrix.metMetaNetXID);
                    excl_mets = {'MNXM1';'MNXM2'};
                    missingidx = ismember(excl_mets,diffmets);
                    missingmets = excl_mets(missingidx);
                    missingidx2 = ~ismember(diffmets,excl_mets);
                    extraidx = ismember(excl_mets,diffmets2);
                    extramets = excl_mets(extraidx);
                    extraidx2 = ~ismember(diffmets2,excl_mets);
                    
                    %Find any missing/extra mets (other than H+/H2O)
                    if any(missingidx2) && ~isempty(diffmets)
                        for l = 1:length(missingidx2)
                            if missingidx2(l) && length(diffmets) > 1
                                arrayidx = find(cellfun('isempty', missingrxnmets),1);
                                missingrxnmets(arrayidx,:) = [model.rxns(i), join([diffmets(l),'is missing from rxn in the model'])];
                            elseif missingidx2(l) && length(diffmets) == 1
                                arrayidx = find(cellfun('isempty', missingrxnmets),1);
                                missingrxnmets(arrayidx,:) = [model.rxns(i), join([diffmets,'is missing from rxn in the model'])];
                            end
                        end
                    elseif any(extraidx2) && ~isempty(diffmets2)
                        for l = 1:length(extraidx2)
                            if extraidx2(l) && length(diffmets2) > 1
                                arrayidx = find(cellfun('isempty', extrarxnmets),1);
                                extrarxnmets(arrayidx,:) = [model.rxns(i), join([diffmets2(l),'is extra in the model rxn'])];
                            elseif extraidx2(l) && length(diffmets2) == 1
                                arrayidx = find(cellfun('isempty', extrarxnmets),1);
                                extrarxnmets(arrayidx,:) = [model.rxns(i), join([diffmets2,'is extra in the model rxn'])];
                            end
                        end
                        
                    else
                        %add/remove H+/H2O into reaction if equation is balanced in MNX database
                        for j = 1:length(excl_mets)
                            if missingidx(j) && isequal(MNXreac_prop(rxnIdx,4),{'TRUE'})
                                mets_MNX(length(mets_MNX) + 1,1) = excl_mets(j);
                                
                                [~,idx] = ismember(mets_MNX,matrix.metMetaNetXID);
                                matrix.metMetaNetXID = matrix.metMetaNetXID(idx);
                                matrix.metcoef = matrix.metcoef(idx);
                                coef_idx = ismember(excl_mets(j),matrix.metMetaNetXID);
                                coef_idx2 = find(ismember(matrix.metMetaNetXID, excl_mets(j)));
                                
                                %Check whether model and MNX database rxn is in the same direction
                                for o = 1:length(mets_coef)
                                    if ~ismember(mets_MNX(o),excl_mets)
                                        if abs(mets_coef(o)) == abs(matrix.metcoef(o)) && mets_coef(o) ~= matrix.metcoef(o)
                                            matrix.metcoef = matrix.metcoef.*(-1);
                                            break
                                        end
                                    end
                                end
                                
                                %check and add coef for H+/H2O
                                if coef_idx && isequal(MNXreac_prop(rxnIdx,4),{'TRUE'})
                                    rxn_comp = unique(mets_comp);
                                    %only change coefficient if there is only 1 compartment
                                    if length(rxn_comp) == 1
                                        mets_locate = find(ismember(model.metMetaNetXID,excl_mets(j)));
                                        metComps2 = modelR.metComps;
                                        mets_comp2 = metComps2(mets_temp);
                                        samemets_comp = find(ismember(modelR.metComps,mets_comp2(1)));
                                        model.S(intersect(mets_locate,samemets_comp),i) = matrix.metcoef(coef_idx2);
                                        coef = model.S(:,i);
                                        mets_temp = (find(model.S(:,i)));
                                        mets_MNX  = model.metMetaNetXID(mets_temp);
                                        mets_coef = coef(mets_temp);
                                        
                                        [~,idx] = ismember(mets_MNX,matrix.metMetaNetXID);
                                        matrix.metMetaNetXID = matrix.metMetaNetXID(idx);
                                        matrix.metcoef = matrix.metcoef(idx);
                                        
                                        %check coef
                                        if isequal(mets_coef,matrix.metcoef) || isequal(mets_coef,matrix.metcoef.*(-1))
                                            arrayidx = find(cellfun('isempty', changeMNXM),1);
                                            changeMNXM(arrayidx,:) = [model.rxns(i),join([excl_mets(j), 'is added with coefficient of',num2str(matrix.metcoef(coef_idx2)),'to balance eqn'])];
                                            model.rxnNotes(i) = join([model.rxnNotes(i),'|',excl_mets(j),'added (PR #222)']);
                                        else
                                            arrayidx = find(cellfun('isempty', coefNoMatch),1);
                                            coefNoMatch(arrayidx,:) = [model.rxns(i),join([excl_mets(j), 'is added with coefficient of',num2str(matrix.metcoef(coef_idx2)),'to balance eqn but coef still does not match'])];
                                        end
                                    else
                                        arrayidx = find(cellfun('isempty', uncertainComp),1);
                                        uncertainComp(arrayidx,:) = [model.rxns(i), join([excl_mets(j),'needs to be added, check which compartment',excl_mets(j),'belongs to in the reaction'])];
                                    end
                                else
                                end
                                
                            elseif extraidx(j) && isequal(MNXreac_prop(rxnIdx,4),{'TRUE'})
                                extramets_comp = mets_comp(find(ismember(mets_MNX,extramets)));
                                extramets_locate = find(ismember(model.metMetaNetXID,extramets));
                                samemets_comp = find(ismember(modelR.metComps,find(ismember(modelR.comps,extramets_comp))));
                                model.S(intersect(extramets_locate,samemets_comp),i) = 0;
                                coef = model.S(:,i);
                                mets_temp = (find(model.S(:,i)));
                                mets_MNX  = model.metMetaNetXID(mets_temp);
                                mets_coef = coef(mets_temp);
                                [~,idx] = ismember(mets_MNX,matrix.metMetaNetXID);
                                matrix.metMetaNetXID = matrix.metMetaNetXID(idx);
                                matrix.metcoef = matrix.metcoef(idx);
                                
                                %check coef
                                if isequal(mets_coef,matrix.metcoef) || isequal(mets_coef,matrix.metcoef.*(-1))
                                    arrayidx = find(cellfun('isempty', changeMNXM),1);
                                    changeMNXM(arrayidx,:) = [model.rxns(i),join([excl_mets(j),'is removed to balance eqn'])];
                                    model.rxnNotes(i) = join([model.rxnNotes(i),'|',excl_mets(j),'is removed to match the balanced equation in MetaNetX database (PR #222)']);
                                else
                                    arrayidx = find(cellfun('isempty', coefNoMatch),1);
                                    coefNoMatch(arrayidx,:) = [model.rxns(i),join([excl_mets(j),'is removed to balance eqn, but model and MNX mets coef still do not match'])];
                                end
                            end
                        end
                    end
                    
                elseif ~any(ismember(mets_MNX,matrix.metMetaNetXID)) 
                    arrayidx = find(cellfun('isempty', allmetsDiff),1);
                    allmetsDiff(arrayidx,:) = [model.rxns(i), join(['All mets in', model.rxns(i), 'of the model are different from database'])];
                    
                elseif all(ismember(mets_MNX,matrix.metMetaNetXID)) && all(ismember(matrix.metMetaNetXID, mets_MNX))
                    %filter transport rxns in the model/rxns with repeated mets
                    if length(mets_MNX) == 2 && ismember(mets_MNX(1),mets_MNX(2))
                    elseif length(unique(mets_MNX)) == 2 && ismember('MNXM1',mets_MNX) %for specific H+ transport reaction
                        %rearrange matrix.metMetaNetXID to match order of mets_MNX
                        %note that there are repeated data in the array
                        [~,idxA] = ismember(mets_MNX,matrix.metMetaNetXID);
                        [~,idxB] = ismember(cellstr(num2str(mets_coef)), cellstr(num2str(matrix.metcoef)));
                        matrix.metMetaNetXID = matrix.metMetaNetXID(idxA);
                        matrix.metcoef = matrix.metcoef(idxB);
                    elseif length(unique(mets_MNX)) < length(mets_MNX) && (~ismember('MNXM1',mets_MNX) || ~ismember('MNXM2',mets_MNX))
                        [~,idxA] = ismember(mets_MNX,matrix.metMetaNetXID);
                        [~,idxB] = ismember(cellstr(num2str(mets_coef)), cellstr(num2str(matrix.metcoef)));
                        matrix.metMetaNetXID = matrix.metMetaNetXID(idxA);
                        matrix.metcoef = matrix.metcoef(idxB);
                    else
                        [~,idx] = ismember(mets_MNX,matrix.metMetaNetXID);
                        matrix.metMetaNetXID = matrix.metMetaNetXID(idx);
                        matrix.metcoef = matrix.metcoef(idx);
                    end
                    
                    if isequal(mets_coef,matrix.metcoef) || isequal(mets_coef,matrix.metcoef.*(-1))
                        arrayidx = find(cellfun('isempty', balanced),1);
                        balanced(arrayidx,:) = [model.rxns(i),'Both mets and coef in model and MNX database matches'];
                    else
                        excl_mets = {'MNXM1';'MNXM2'};
                        %check whether model and MNX database rxn is in the same direction
                        for o = 1:length(mets_coef)
                            if ~ismember(mets_MNX(o),excl_mets)
                                if abs(mets_coef(o)) == abs(matrix.metcoef(o)) && mets_coef(o) ~= matrix.metcoef(o)
                                    matrix.metcoef = matrix.metcoef.*(-1);
                                    break
                                end
                            end
                        end
                        for k = 1:length(excl_mets)
                            coef_idx = ismember(excl_mets(k),matrix.metMetaNetXID);
                            coef_idx2 = find(ismember(matrix.metMetaNetXID, excl_mets(k)));
                            %check and correct unbalanced H+/H2O
                            if coef_idx && isequal(MNXreac_prop(rxnIdx,4),{'TRUE'})
                                for m = 1:length(coef_idx2)
                                    if abs(matrix.metcoef(coef_idx2(m))) > abs(mets_coef(coef_idx2(m))) && mets_coef(coef_idx2(m)) < 0
                                        model.rxnNotes(i) = join([model.rxnNotes(i),'| Coefficient of',excl_mets(k),'modified from',num2str(mets_coef(coef_idx2(m))),'to',num2str(matrix.metcoef(coef_idx2(m))),'(PR #222)']);
                                        mets_coef(coef_idx2(m)) = matrix.metcoef(coef_idx2(m));
                                        mets_locate = find(ismember(model.metMetaNetXID,excl_mets(k)));
                                        metComps2 = modelR.metComps;
                                        mets_comp2 = metComps2(mets_temp);
                                        samemets_comp = find(ismember(modelR.metComps,mets_comp2(coef_idx2(m))));
                                        model.S(intersect(mets_locate,samemets_comp),i) = mets_coef(coef_idx2(m));
                                        
                                    elseif abs(matrix.metcoef(coef_idx2(m))) > abs(mets_coef(coef_idx2(m))) &&  mets_coef(coef_idx2(m)) > 0
                                        model.rxnNotes(i) = join([model.rxnNotes(i),'| Coefficient of',excl_mets(k),'modified from',num2str(mets_coef(coef_idx2(m))),'to',num2str(matrix.metcoef(coef_idx2(m))),'(PR #222)']);
                                        mets_coef(coef_idx2(m)) = matrix.metcoef(coef_idx2(m));
                                        mets_locate = find(ismember(model.metMetaNetXID,excl_mets(k)));
                                        metComps2 = modelR.metComps;
                                        mets_comp2 = metComps2(mets_temp);
                                        samemets_comp = find(ismember(modelR.metComps,mets_comp2(coef_idx2(m))));
                                        model.S(intersect(mets_locate,samemets_comp),i) = mets_coef(coef_idx2(m));
                                        
                                    elseif abs(matrix.metcoef(coef_idx2(m))) < abs(mets_coef(coef_idx2(m))) && mets_coef(coef_idx2(m)) < 0
                                        model.rxnNotes(i) = join([model.rxnNotes(i),'| Coefficient of',excl_mets(k),'modified from',num2str(mets_coef(coef_idx2(m))),'to',num2str(matrix.metcoef(coef_idx2(m))),'(PR #222)']);
                                        mets_coef(coef_idx2(m)) = matrix.metcoef(coef_idx2(m));
                                        mets_locate = find(ismember(model.metMetaNetXID,excl_mets(k)));
                                        metComps2 = modelR.metComps;
                                        mets_comp2 = metComps2(mets_temp);
                                        samemets_comp = find(ismember(modelR.metComps,mets_comp2(coef_idx2(m))));
                                        model.S(intersect(mets_locate,samemets_comp),i) = mets_coef(coef_idx2(m));
                                        
                                    elseif abs(matrix.metcoef(coef_idx2(m))) < abs(mets_coef(coef_idx2(m))) && mets_coef(coef_idx2(m)) > 0
                                        model.rxnNotes(i) = join([model.rxnNotes(i),'| Coefficient of',excl_mets(k),'modified from',num2str(mets_coef(coef_idx2(m))),'to',num2str(matrix.metcoef(coef_idx2(m))),'(PR #222)']);
                                        mets_coef(coef_idx2(m)) = matrix.metcoef(coef_idx2(m));
                                        mets_locate = find(ismember(model.metMetaNetXID,excl_mets(k)));
                                        metComps2 = modelR.metComps;
                                        mets_comp2 = metComps2(mets_temp);
                                        samemets_comp = find(ismember(modelR.metComps,mets_comp2(coef_idx2(m))));
                                        model.S(intersect(mets_locate,samemets_comp),i) = mets_coef(coef_idx2(m));
                                        
                                    elseif abs(matrix.metcoef(coef_idx2(m))) == abs(mets_coef(coef_idx2(m))) && matrix.metcoef(coef_idx2(m)) ~= mets_coef(coef_idx2(m))
                                        model.rxnNotes(i) = join([model.rxnNotes(i),'| Coefficient of',excl_mets(k),'modified from',num2str(mets_coef(coef_idx2(m))),'to',num2str(matrix.metcoef(coef_idx2(m))),'(PR #222)']);
                                        mets_coef(coef_idx2(m)) = matrix.metcoef(coef_idx2(m));
                                        mets_locate = find(ismember(model.metMetaNetXID,excl_mets(k)));
                                        metComps2 = modelR.metComps;
                                        mets_comp2 = metComps2(mets_temp);
                                        samemets_comp = find(ismember(modelR.metComps,mets_comp2(coef_idx2(m))));
                                        model.S(intersect(mets_locate,samemets_comp),i) = mets_coef(coef_idx2(m));
                                    else
                                    end
                                end
                            end
                        end
                        if isequal(mets_coef,matrix.metcoef) || isequal(mets_coef,matrix.metcoef.*(-1))
                            arrayidx = find(cellfun('isempty', changeMNXM),1);
                            changeMNXM(arrayidx,:) = [model.rxns(i),'Coef for H+/H2O has been modified (as indicated in model.rxnNotes) to balance rxn in model'];
                        else
                            arrayidx = find(cellfun('isempty', coefNoMatch),1);
                            coefNoMatch(arrayidx,:) = [model.rxns(i),'Mets in model and MNX database matches, but coef does not match'];
                        end
                    end
                end
            end
        end
    end
    
    %remove empty cells from cell arrays
    arrayidx = ~cellfun('isempty',balanced(:,1));
    strtemp = balanced(:,1);
    strtemp2 = balanced(:,2);
    balanced = [strtemp(arrayidx),strtemp2(arrayidx)];
    
    arrayidx = ~cellfun('isempty',changeMNXM(:,1));
    strtemp = changeMNXM(:,1);
    strtemp2 = changeMNXM(:,2);
    changeMNXM = [strtemp(arrayidx),strtemp2(arrayidx)];
    
    arrayidx = ~cellfun('isempty',noMNXRID(:,1));
    strtemp = noMNXRID(:,1);
    strtemp2 = noMNXRID(:,2);
    noMNXRID = [strtemp(arrayidx),strtemp2(arrayidx)];
    
    arrayidx = ~cellfun('isempty',noCoef(:,1));
    strtemp = noCoef(:,1);
    strtemp2 = noCoef(:,2);
    noCoef = [strtemp(arrayidx),strtemp2(arrayidx)];
    
    arrayidx = ~cellfun('isempty',allmetsDiff(:,1));
    strtemp = allmetsDiff(:,1);
    strtemp2 = allmetsDiff(:,2);
    allmetsDiff = [strtemp(arrayidx),strtemp2(arrayidx)];
    
    arrayidx = ~cellfun('isempty',coefNoMatch(:,1));
    strtemp = coefNoMatch(:,1);
    strtemp2 = coefNoMatch(:,2);
    coefNoMatch = [strtemp(arrayidx),strtemp2(arrayidx)];
    
    arrayidx = ~cellfun('isempty',missingrxnmets(:,1));
    strtemp = missingrxnmets(:,1);
    strtemp2 = missingrxnmets(:,2);
    missingrxnmets = [strtemp(arrayidx),strtemp2(arrayidx)];
    
    arrayidx = ~cellfun('isempty',extrarxnmets(:,1));
    strtemp = extrarxnmets(:,1);
    strtemp2 = extrarxnmets(:,2);
    extrarxnmets = [strtemp(arrayidx),strtemp2(arrayidx)];
    
    %clear all variables except outputs of script
    clearvars -except model noMNXRID noCoef coefNoMatch...
        missingrxnmets balanced changeMNXM
    
    %remove reac_prop.tsv from current directory
    delete('reac_prop.tsv');
    
    % Remove leading '| ' in notes that were previously empty
    model.rxnNotes = regexprep(model.rxnNotes,'^\| ','');
    model.metNotes = regexprep(model.metNotes,'^\| ','');
else
    warning('reac_prop.tsv is not downloaded into current directory, S matrix is not checked');
end

%save model
saveYeastModel(model);

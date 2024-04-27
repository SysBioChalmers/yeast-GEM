function model = addYMDBconcentrations(model,downloadYMDB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = addYMDBconcentrations(model)
% Loads the YMDB metabolite concentrations as available from
% data/databases/YMDBconcData.csv, and assign them in a model.metConcs
% structure, which contains three vectors with 'mean', 'max' and 'min'
% concentrations for each metabolite. Data is matched by KEGG identifier,
% CHEBI and metabolite name. Note that the .metConcs fields are not kept
% when the model is exported in any model format.
%
% If downloadYMDB is set as true, the data will be downloaded from YMDB and
% processed, otherwise the existing file is used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2
    downloadYMDB = false;
end
if downloadYMDB
    downloadYMDBdata();
elseif ~exist('../../data/databases/YMDBconcentrations.csv','file')
    downloadYMDBdata();
end

% Load YMDB data
data = readtable('../../data/databases/YMDBconcentrations.csv');
data = table2cell(data);
metConcs.mean = nan(length(model.metNames), 1);
metConcs.max = nan(length(model.metNames), 1);
metConcs.min = nan(length(model.metNames), 1);

met.name = model.metNames;
met.miriam = extractMiriam(model.metMiriams,{'kegg.compound','chebi'});
met.kegg = met.miriam(:,1);
met.kegg(cellfun(@isempty,met.kegg)) = {'###'}; % Prevent matching of empty identifiers
met.chebi = met.miriam(:,2);
met.chebi(cellfun(@isempty,met.chebi)) = {'###'};

[lia,locb]=ismember(met.chebi,data(:,3));
metConcs.mean(lia,1) = data{locb(lia),5};
metConcs.max(lia,1) = data{locb(lia),6};
metConcs.min(lia,1) = data{locb(lia),7};

[lia,locb]=ismember(met.kegg,data(:,4));
metConcs.mean(lia,1) = data{locb(lia),5};
metConcs.max(lia,1) = data{locb(lia),6};
metConcs.min(lia,1) = data{locb(lia),7};

[lia,locb]=ismember(lower(met.name),lower(data(:,2)));
metConcs.mean(lia,1) = data{locb(lia),5};
metConcs.max(lia,1) = data{locb(lia),6};
metConcs.min(lia,1) = data{locb(lia),7};

fprintf('In total, %d metabolites have concentration data.\n', length(find(~isnan(metConcs.mean))));
model.metConcs = metConcs;
end

function downloadYMDBdata()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% downloadYMDBdata
% Downloads metabolite concentrations from YMDB and stores all data in
% data/databases/YMDBconcData.csv. Uses ProgressBar2 from RAVEN for
% reporting progress.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

raw = transpose(char(webread('https://www.ymdb.ca/compounds.csv?compound_quantified=1')));
IDs = strsplit(raw,'\n')';
IDs(1) = [];
IDs(cellfun(@isempty,IDs))=[];
IDs = regexp(IDs,'^YMDB\d{5}','match');

ymdbNum = numel(IDs);
YMDBmax = zeros(ymdbNum,1);
YMDBmin = zeros(ymdbNum,1);
YMDBavg = zeros(ymdbNum,1);
YMDBdat = cell(ymdbNum,1);
metName = cell(ymdbNum,1);
KEGGid = cell(ymdbNum,1);
KEGGid(:) = {''};
CHEBIid = cell(ymdbNum,1);
CHEBIid(:) = {''};
showProgress = exist('ProgressBar2','file');
if showProgress
    PB = ProgressBar2(ymdbNum,'Downloading from YMDB','gui');
end
parfor i=1:ymdbNum
    cmp = IDs{i};
    url = ['https://www.ymdb.ca/compounds/',cmp{1}];
    success = false;
    raw = ''; out='';
    try
        raw = webread(url);
        success = true;
    catch
        fprintf('Failed to download %s\n',IDs{i}{1})
    end
    if success
        concStart = strfind(raw,'Intracellular Concentration');
        concEnd   = strfind(raw,'Extracellular Concentration');
        if numel(concStart)>1
            out=raw(concStart(2):concEnd);
            out=strsplit(out,'<td>');
            numData = (numel(out)-2)/5;
        else
            numData = -1000;
        end
        if numData < 1 && numData > -1000
            concEnd   = strfind(raw,'Conversion Details Here');
            out=raw(concStart(2):concEnd);
            out=strsplit(out,'<td>');
            numData = (numel(out)-2)/5;
        end
        if numData > 0
            out(1)=[];
            concStart = 1:5:((numData-1)*5+1);
            concData = out(concStart)';
            concData = regexprep(concData,' <\/td>','');
            measData = regexp(concData,'^[\d\.]+','match');
            measData = cellstr([measData{:}]);
            measData = cellfun(@str2num,measData);
            errorData = regexp(concData,' [\d\.]+','match');
            errorData = cellstr([errorData{:}]);
            errorData = cellfun(@str2num,errorData);
            maxConc = max(measData+errorData);
            minConc = min(measData-errorData);
            meanConc = mean(measData);

            concData = regexprep(concData,'&plusmn;','+/-');
            concData = regexprep(concData,'&#181;','u');
            YMDBdat{i} = strjoin(concData,'; ');

            metName(i) = regexp(raw,'Name<\/th><td>(.+)<\/td><\/tr><tr><th>Species','tokens');

            keggQ = '.*http\:\/\/www\.genome\.jp\/dbget-bin\/www_bget\?cpd\:(C\d{5}).*';
            tmpK = regexp(raw,keggQ,'tokens');
            if ~isempty(tmpK)
                KEGGid(i,1) = tmpK{1};
            end

            tmpC = regexp(raw,'chebiId=(\d+)','tokens');
            if ~isempty(tmpC)
                CHEBIid(i,1) = append('CHEBI:',tmpC{1});
            end

            YMDBmax(i) = maxConc;
            YMDBmin(i) = minConc;
            YMDBavg(i) = meanConc;
        else
            fprintf('Failed to extract %s\n',IDs{i}{1})
        end
    end
    if showProgress
        count(PB)
    end
end

noData = cellfun(@isempty,metName);

metConcs = [IDs(~noData),metName(~noData),CHEBIid(~noData),KEGGid(~noData),num2cell(YMDBavg(~noData)),num2cell(YMDBmax(~noData)),num2cell(YMDBmin(~noData)),YMDBdat(~noData)];
metConcs = [{'IDs','metNames','chebi','kegg','mean','max','min','concs'};metConcs];
metConcs = cell2table(metConcs);
writetable(metConcs, '../../data/databases/YMDBconcentrations.csv', 'FileType', 'text', 'WriteVariableNames',false);
end

function increaseVersion(bumpType)
% increaseVersion
%   Upgrades the model to a new version. Run this function after merging
%   changes to the main branch for making a new release.
%
%   bumpType    One of the following 3 strings: "major", "minor" or
%               "patch", indicating the type of increase of version to be
%               performed.
%
%   Usage: increaseVersion(bumpType)
%

%Check if in main:
[~,currentBranch] = system('git rev-parse --abbrev-ref HEAD');
if ~strcmp(strtrim(currentBranch),'main')
    error('ERROR: not in main')
end

%Bump version number:
fid = fopen('../version.txt','r');
oldVersion = fgetl(fid);
fclose(fid);
oldVersion = str2double(strsplit(oldVersion,'.'));
newVersion = oldVersion;
switch bumpType
    case 'major'
        newVersion(1) = newVersion(1) + 1;
        newVersion(2) = 0;
        newVersion(3) = 0;
    case 'minor'
        newVersion(2) = newVersion(2) + 1;
        newVersion(3) = 0;
    case 'patch'
        newVersion(3) = newVersion(3) + 1;
    otherwise
        error('ERROR: invalid input. Use "major", "minor" or "patch"')
end
newVersion = num2str(newVersion,'%d.%d.%d');

%Check if history has been updated:
fid     = fopen('../history.md','r');
history = fscanf(fid,'%s');
fclose(fid);
if ~contains(history,['yeast' newVersion ':'])
    error('ERROR: update history.md first')
end

%Load model:
disp('Loading model file')
model = readYAMLmodel('../model/yeast-GEM.yml');

%Run tests
cd modelTests
disp('Running gene essentiality analysis')
[new.accuracy,new.tp,new.tn,new.fn,new.fp] = essentialGenes(model,true);
disp('Run growth analysis')
new.R2=growth(model,true);
disp('Run anaerobic flux analysis')
cd ../otherChanges
modelAnaerobic = anaerobicModel(model);
cd ../modelTests
new.R2anaerobic = anaerobic_flux_predictions(modelAnaerobic);
disp('Run anaerobic exchange rate analysis')
[new.exchange, new.exchangeTable] = plotAnaerobic(modelAnaerobic,false);

cd ..
copyfile('../README.md','backup.md')
fin  = fopen('backup.md','r');
fout = fopen('../README.md','w');
searchStats1 = '^(- Accuracy\: )0\.\d+';
searchStats2 = '^(- True non-essential genes\: )\d+';
searchStats3 = '^(- True essential genes\: )\d+';
searchStats4 = '^(- False non-essential genes\: )\d+';
searchStats5 = '^(- False essential genes\: )\d+';
newStats1 = ['$1' num2str(new.accuracy,'%.3f')];
newStats2 = ['$1' num2str(numel(new.tp))];
newStats3 = ['$1' num2str(numel(new.tn))];
newStats4 = ['$1' num2str(numel(new.fp))];
newStats5 = ['$1' num2str(numel(new.fn))];

searchStats6 = '^(- Correlation coefficient R<sup>2<\/sup>\: )0\.\d+';
newStats6 = ['$1' num2str(new.R2,'%.3f')];

while ~feof(fin)
    str = fgets(fin);
    inline = regexprep(str,searchStats1,newStats1);
    inline = regexprep(inline,searchStats2,newStats2);
    inline = regexprep(inline,searchStats3,newStats3);
    inline = regexprep(inline,searchStats4,newStats4);
    inline = regexprep(inline,searchStats5,newStats5);
    inline = regexprep(inline,searchStats6,newStats6);
    inline = unicode2native(inline,'UTF-8');
    fwrite(fout,inline);
end
fclose('all');
delete('backup.md');

%Update the test results summary. Only the two metadata lines and the
%numbers in the table are rewritten, so the surrounding text describing
%what each test does can be edited freely. The Python parity gate reads
%this same file, so it must not be allowed to go stale.
resultsFile = '../data/testResults/README.md';
copyfile(resultsFile,'backup.md')
fin  = fopen('backup.md','r');
fout = fopen(resultsFile,'w');
searchRes = {'^(- Model version\: )[^\n]*', ...
             '^(- Software\: )[^\n]*', ...
             '^(\| Growth prediction R2 \| )[\d.eE+-]+', ...
             '^(\| Anaerobic flux prediction R2 \| )[\d.eE+-]+', ...
             '^(\| Anaerobic exchange mean relative error \| )[\d.eE+-]+', ...
             '^(\| Anaerobic exchange within error \| )[\d.eE+-]+', ...
             '^(\| Ammonium per ATPase \| )[\d.eE+-]+', ...
             '^(\| Gene essentiality accuracy \| )[\d.eE+-]+', ...
             '^(\| True non-essential genes \| )\d+', ...
             '^(\| True essential genes \| )\d+', ...
             '^(\| False non-essential genes \| )\d+', ...
             '^(\| False essential genes \| )\d+'};
newRes    = {['$1' newVersion], ...
             ['$1MATLAB ' version ', RAVEN ' getToolboxVersion('RAVEN','ravenCobraWrapper.m')], ...
             ['$1' num2str(new.R2,'%.16g')], ...
             ['$1' num2str(new.R2anaerobic,'%.16g')], ...
             ['$1' num2str(new.exchange.meanRelativeError,'%.16g')], ...
             ['$1' num2str(new.exchange.fractionWithinError,'%.16g')], ...
             ['$1' num2str(new.exchange.ammoniumPerATPase,'%.16g')], ...
             ['$1' num2str(new.accuracy,'%.16g')], ...
             ['$1' num2str(numel(new.tp))], ...
             ['$1' num2str(numel(new.tn))], ...
             ['$1' num2str(numel(new.fp))], ...
             ['$1' num2str(numel(new.fn))]};

%One row per measured exchange rate. The measured value and its error come
%from the data file and do not change, but the whole row is rewritten so
%that the row and the prediction cannot fall out of step.
exLabels = regexprep(new.exchangeTable.rxnName,' (exchange|pseudoreaction)$','');
for i = 1:numel(exLabels)
    if new.exchangeTable.withinError(i)
        within = 'yes';
    else
        within = 'no';
    end
    searchRes{end+1} = ['^\| ' exLabels{i} ' \|[^\n]*']; %#ok<AGROW>
    newRes{end+1} = sprintf('| %s | %g +/- %g | %.4f | %s |', ...
        exLabels{i}, new.exchangeTable.measured(i), ...
        new.exchangeTable.stdev(i), new.exchangeTable.predicted(i), within); %#ok<AGROW>
end

while ~feof(fin)
    inline = fgets(fin);
    for i = 1:numel(searchRes)
        inline = regexprep(inline,searchRes{i},newRes{i});
    end
    inline = unicode2native(inline,'UTF-8');
    fwrite(fout,inline);
end
fclose('all');
delete('backup.md');

%Include tag and save model:
disp('Write model files')
model.id = ['yeastGEM_v' newVersion];
model.version = newVersion;
saveYeastModel(model,true,true,true)   %only save if model can grow

%Stage the binary model files. They are listed in .gitignore, so that they
%cannot reach develop by accident, and are therefore force-added here on
%main where they do belong. This replaces stripping the two entries out of
%.gitignore: that left main and develop permanently disagreeing about those
%lines, which made every merge between the two branches conflict.
system('git add -f ../model/yeast-GEM.mat ../model/yeast-GEM.xlsx');

%Check for any unexpected file changes
[~,diff]   = system('git diff --numstat');
diff   = strsplit(diff,'\n');
change = false;
for i = 1:length(diff)
    diff_i = strsplit(diff{i},'\t');
    if length(diff_i) == 3
        switch diff_i{3}
            case 'model/yeast-GEM.xml'
                %.xml file: 4 lines should be added & 4 lines should be
                %deleted (2 with version information, 2 with current date)
                if eval([diff_i{1} ' > 4']) || eval([diff_i{2} ' > 4'])
                    disp(['NOTE: File ' diff_i{3} ' is changing more than expected'])
                    change = true;
                end
            case 'model/yeast-GEM.yml'
                %.yml file: 3 lines should be added & 3 lines should be deleted
                %(2 with version information, 1 with current date)
                if eval([diff_i{1} ' > 3']) || eval([diff_i{2} ' > 3'])
                    disp(['NOTE: File ' diff_i{3} ' is changing more than expected'])
                    change = true;
                end                
            %The test results are regenerated above, so they are expected
            %to change whenever the model does
            case {'history.md','README.md','model/yeast-GEM.mat','model/yeast-GEM.xlsx', ...
                  'data/testResults/README.md','data/testResults/growth.md', ...
                  'data/testResults/growth.png','data/testResults/essentialGenes.md'}
            otherwise
                disp(['NOTE: File ' diff_i{3} ' is changing'])
                change = true;                
        end
    end
end
if change == true
    error(['Some files are changing from develop. To fix, first update develop, ' ...
        'then merge to main, and try again.'])
end

%Update version file:
fid = fopen('../version.txt','wt');
fprintf(fid,newVersion);
fclose(fid);
end

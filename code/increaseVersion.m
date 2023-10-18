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
[new.accuracy,new.tp,new.tn,new.fn,new.fp] = essentialGenes(model);
disp('Run growth analysis')
new.R2=growth(model);

saveas(gcf,'../../growth.png');

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

%Allow .mat & .xlsx storage:
copyfile('../.gitignore','backup')
fin  = fopen('backup','r');
fout = fopen('../.gitignore','w');
still_reading = true;
while still_reading
  inline = fgets(fin);
  if ~ischar(inline)
      still_reading = false;
  elseif ~startsWith(inline,'*.mat') && ~startsWith(inline,'*.xls')
      fwrite(fout,inline);
  end
end
fclose('all');
delete('backup');

%Include tag and save model:
disp('Write model files')
model.id = ['yeastGEM_v' newVersion];
model.version = newVersion;
saveYeastModel(model,true,true,true)   %only save if model can grow

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
            case {'history.md','README.md','growth.png','model/yeast-GEM.mat'}
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

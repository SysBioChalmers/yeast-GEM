function model = getEarlierModelVersion(version,verbose)
% getEarlierModelVersion
%   Loads an earlier model version from the local git repository history,
%   and removes any mention of its version from model.id.
%   Requires that the yeast-GEM repository is cloned via git.
%
%   version     either 'main' or 'develop' for the latest model file from
%               the specified branch. Can alternatively be a specific model
%               version (e.g. '8.5.1'), to load a the model file from a
%               specific release version. 
%   verbose     if true, the model version (as obtained from model.id) is
%               printed (opt, default true).
%
% Usage:
%   model = getEarlierModelVersion(version,verbose)

if nargin<2
    verbose=true;
end
%Test that git is working
[~,out]=system('git --version');
if ~startsWith(out,'git version ')
    error('git does not seem to be installed')
end

if strcmp(version,'main')
    system('git show main:model/yeast-GEM.xml > _earlierModel.xml');
elseif strcmp(version,'develop')
    system('git show develop:model/yeast-GEM.xml > _earlierModel.xml');
elseif ~isempty(regexp(version,'^\d+\.\d+\.\d+$'))
    tagpath= ['refs/tags/v' version ':model/yeast-GEM.xml'];
    [out,~]=system(['git show ' tagpath ' > _earlierModel.xml']);
    if out==128 % File not found, try older folder structure
        delete '_earlierModel.xml'
        tagpath= ['refs/tags/v' version ':ModelFiles/xml/yeastGEM.xml'];
        [out,~]=system(['git show ' tagpath ' > _earlierModel.xml']);
        if out==128
            delete '_earlierModel.xml'
            error('Cannot find the specified model version')
        end
    end
else
    error('Unclear which earlier model version should be loaded')
end

model=importModel('_earlierModel.xml');
delete '_earlierModel.xml'
if verbose
    disp(['Loaded model version: ''' regexprep(model.id,'yeastGEM_v?','') ''''])
end
end

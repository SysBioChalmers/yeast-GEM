function root = yeastRoot()
% yeastRoot  The yeast-GEM checkout that the analyses use: the folder in the
%   environment variable YEAST_GEM_ROOT, or else the repository that contains
%   this folder.

root = getenv('YEAST_GEM_ROOT');
if isempty(root)
    root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
end
end

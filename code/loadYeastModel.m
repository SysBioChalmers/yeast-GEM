function model = loadYeastModel(filename)
% loadYeastModel
%   Load the yeast-GEM in a MATLAB environment, requires the RAVEN Toolbox.
%   If RAVEN Toolbox is not installed, an attempt is made to load the model
%   via COBRA Toolbox, accompanied with a warning as RAVEN is recommended.
%
% Input:
%   filename    by default, the model is loaded from its location at
%               yeast-GEM/model/yeast-GEM.yml. Alternative model files can
%               be loaded if provided here. (opt, default empty).
%
% Output:
%   model       the yeast-GEM model structure. Reaction, metabolite and
%               gene cross-reference annotation (KEGG, BiGG, ChEBI,
%               MetaNetX, EC codes, UniProt) is stored separately in
%               model/reactions.tsv, model/metabolites.tsv and
%               model/genes.tsv (yeast-GEM#379) and merged back in here
%               via annotateGEM, so the returned model looks the same as
%               before that split either way.
%
%   Usage: model = loadYeastModel(filename)

funcDir = dbstack('-completenames');
funcDir = regexprep(funcDir(1).file,[funcDir(1).name '\.m'],'');

if nargin<1 || isempty(filename)
    filename = fullfile(funcDir,'..','model','yeast-GEM.yml');
end

if endsWith(filename,{'.yml','.yaml'})
    model = readYAMLmodel(filename);
    % The yml only carries this by itself for the default model path -- a
    % filename argument may point at an .xml-derived .yml or some other
    % file that already has full annotation, or at a model/ folder
    % without the tsvs, in which case there is nothing to merge.
    annPath = fullfile(funcDir,'..','model');
    if isequal(filename, fullfile(funcDir,'..','model','yeast-GEM.yml')) ...
            && isfile(fullfile(annPath,'reactions.tsv'))
        model = annotateGEM(model, annPath);
    end
else
    if ~(exist('ravenCobraWrapper.m','file')==2)
        if exist('readCbModel.m','file')==2
            warning(['RAVEN cannot be found. Attempt to load yeast-GEM in '...
                'COBRA format.\n\nNote that it is recommended to have RAVEN '...
                'installed, especially when curating yeast-GEM (see README.md for '...
                'more info).%s'],'')
            model = readCbModel(filename);
        else
            error(['RAVEN cannot be found. See README.md for installation '...
                'instructions.'])
        end
    else
        model = importModel(filename);
    end
    currentDir = pwd;
    cd(fullfile(funcDir,'missingFields'));
    model = loadDeltaG(model);
    cd(currentDir)
end
end

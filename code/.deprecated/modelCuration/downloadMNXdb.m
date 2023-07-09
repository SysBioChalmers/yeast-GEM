function downloadMNXdb(files,mnxPath)
%downloadMNXdb  Download MetaNetX database files if they cannot be found.
%
%   files       string or cell array of strings, with the MetaNetX files to
%               be downloaded. Options: 'chem_xref', 'chem_prop',
%               'reac_xref', 'reac_prop' or 'all'.
%   mnxPath     string of path where MetaNetX reference files are
%               stored. (opt, default to RAVENdir/external/metanetx) To
%               download to current folder, specify pwd().
%
% Usage: downloadMNXdb(files,mnxPath)
%

if nargin<2
    [ST, I]=dbstack('-completenames');
    mnxPath=fileparts(fileparts(fileparts(ST(I).file)));
    mnxPath=fullfile(mnxPath,'external','metanetx');
end

if ischar(files)
    if strcmp(files,'all')
        files={'chem_xref','chem_prop','reac_prop','reac_xref'};
    else
        files={files};
    end
end

for i=1:length(files)
    if ~exist(fullfile(mnxPath,[files{i},'.tsv']), 'file')
        fprintf('File %s.tsv cannot be found and will be downloaded from MetaNetX.org.\n',files{i});
        websave(fullfile(mnxPath,[files{i},'.tsv']),...
            ['https://www.metanetx.org/cgi-bin/mnxget/mnxref/',files{i},'.tsv']);
    end
end
end


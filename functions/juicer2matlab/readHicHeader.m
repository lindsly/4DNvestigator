function [hicHeader] = readHicHeader(fn)
%readHicHeader reads the .hic header information
%
%   Input
%   fn:         full path to .hic file
%
%   Output
%   hicHeader: .hic header information
%
%   Scott Ronquist, scotronq@umich.edu. 1/22/19

%% format juicer tools input info
% find 4DNvestigator juicer tools path
juicerJarDir = mfilename('fullpath');

folderSlash = '\';
if isunix
    folderSlash = '/';
end

juicerJarDirLevels = strfind(juicerJarDir,folderSlash);
juicerJarDir = [juicerJarDir(1:juicerJarDirLevels(end)),'juicer_tools.jar'];

%% get and format .hic header info
% get "read_hic_header.py" file location
juicerReadHic = [juicerJarDir(1:juicerJarDirLevels(end)),...
    '\straw-master\python\read_hic_header.py'];

% run "read_hic_header.py"
[status,cmdout] = system(sprintf('python %s %s',juicerReadHic,fn));

% report an error if status~=0
if status~=0
    error(cmdout)
end

% parse header info
[hicHeader] = parseHicHeaderInfo(cmdout);

end


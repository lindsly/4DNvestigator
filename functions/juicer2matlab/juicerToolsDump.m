function [juicerOut,hicHeader] = juicerToolsDump(norm3d,norm1d,fn,chr1,chr2,bpFrag,binSize,fnOut,headerFlag)
%juicer_tools_dump_mat MATLAB version of juicertools "dump" to read .hic
%   juicer_tools_dump_mat reads in .hic files and outputs contact matrices
%   for the input parameters. Functions inputs follow juicer_tools dump
%   command input.
%
%   Reference:
%   https://github.com/theaidenlab/juicer/wiki/Data-Extraction
%
%   Input
%   norm1d: 2D normalization [observed/oe]
%   norm1d: 1D normalization [NONE/VC/VC_SQRT/KR]
%   fn:     hicFile(s) location
%   chr1:   chromosome # (eg 1-22,X,Y in human)
%   chr2:   chromosome # (eg 1-22,X,Y in human)
%   bpFrag: bin units [BP/FRAG] (FRAG dependant on RE)
%   binSize: bin size (ie 1E5 for 100kb resolution)
%   fnOut:  temporary file name to output. Reccomended that this not input
%           this, and let it be default (this will create a temporary .txt file in 
%           your wd, which will be deleted automatically)
%   headerFlag: flag for whether or not to output hic header (requires python)
%
%   Output
%   juicerOut: juicer_tools dump output (N x 3 matrix)
%   hicHeader: MATLAB structure with .hic header information
%
%   example:
%
%   Scott Ronquist, scotronq@umich.edu. 6/27/18, 12/12/18

if ~exist('fnOut','var')||isempty(fnOut);fnOut = 'juicer_temp.txt';end
if ~exist('headerFlag','var')||isempty(headerFlag);headerFlag = 0;end

%% num2str if necessary
if isnumeric(chr1);chr1=num2str(chr1);end
if isnumeric(chr2);chr2=num2str(chr2);end
if isnumeric(binSize);binSize=num2str(binSize);end

%% format juicer tools input info
% find GSFAT juicer tools path
juicerJarDir = mfilename('fullpath');

folderSlash = '\';
if isunix
    folderSlash = '/';
end

juicerJarDirLevels = strfind(juicerJarDir,folderSlash);
juicerJarDir = [juicerJarDir(1:juicerJarDirLevels(end)),'juicer_tools.jar'];

%% run juicer Dump
[status,cmdout] = system(sprintf('java -jar "%s" dump %s %s %s %s %s %s %s %s',...
    juicerJarDir,norm3d,norm1d,fn,chr1,chr2,bpFrag,binSize,fnOut));

% report an error if status~=0
if status~=0
    error(cmdout)
end

%% get .hic file header
% only extract if specified by function input
if headerFlag
    % get "read_hic_header.py" file location
    juicerReadHic = [juicerJarDir(1:juicerJarDirLevels(end)),...
        '\straw-master\python\read_hic_header.py'];
    
    % run "read_hic_header.py"
    [status2,cmdout2] = system(sprintf('python %s %s',juicerReadHic,fn));
    
    % report an error if status~=0
    if status2~=0
        error(cmdout2)
    end
    
    % parse header info
    [hicHeader] = parseHicHeaderInfo(cmdout2);
else
    hicHeader = struct;
end

%% format hic dump matrix data
if ~isempty(fnOut)
    temp = readtable(fnOut,'Delimiter','\t');
    juicerOut = temp{:,:};
    if isempty(juicerOut)
        juicerOut = zeros(1,3);
    elseif ~strcmpi(chr1,'ALL')
        juicerOut(:,1:2) = round(juicerOut(:,1:2)/str2num(binSize));
    end
    if strcmp(fnOut,'juicer_temp.txt')
        delete juicer_temp.txt
    end
else
    lineLocs = find(ismember(cmdout, char([10 13])));
    lineLocs = [0,lineLocs(1:end-1)];
    juicerOut = zeros(length(lineLocs),3);
    for i = 1:length(lineLocs)-1
        fprintf('formatting juicer: %.4f\n',i/(length(lineLocs)-1))
        juicerOut(i,:) = str2num(cmdout(lineLocs(i)+1:lineLocs(i+1)-1));
    end
    juicerOut(:,1:2) = round(juicerOut(:,1:2)/str2num(binSize));
end

end


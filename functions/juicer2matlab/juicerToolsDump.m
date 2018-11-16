function [juicerOut] = juicerToolsDump(norm2d,norm1d,fn,chr1,chr2,bpFrag,binSize,fnOut)
%juicer_tools_dump_mat MATLAB version of juicertools "dump" to read .hic
%   juicer_tools_dump_mat reads in .hic files and outputs contact matrices
%   for the input parameters. Functions inputs follow juicer_tools dump
%   command input.
%
%   Reference:
%   https://github.com/theaidenlab/juicer/wiki/Data-Extraction
%
%   Inputs:
%   norm_2d: 2D normalization [observed/oe]
%   norm_1d: 1D normalization [NONE/VC/VC_SQRT/KR]
%   fn: hicFile(s) location
%   chr1: chromosome # (eg 1-22,X,Y in human)
%   chr2: chromosome # (eg 1-22,X,Y in human)
%   bp_frag: bin units [BP/FRAG] (FRAG dependant on RE)
%   bin_size: bin size (ie 1E5 for 100kb resolution)
%   fn_out: temporary file name to output. Reccomended that this not input
%   this, and let it be default (this will create a temporary .txt file in 
%   your wd, which will be deleted automatically)
%
%   juicer_out: juicer_tools dump output (N x 3 matrix)
%
%   example:
%
%   Scott Ronquist, 6/27/18

if ~exist('fnOut','var')||isempty(intraFfnOutlag);fnOut = 'juicer_temp.txt';end

%% num2str if necessary
if isnumeric(chr1);chr1=num2str(chr1);end
if isnumeric(chr2);chr2=num2str(chr2);end
if isnumeric(binSize);binSize=num2str(binSize);end

%% call juicer tools
% find GSFAT juicer tools path
juicerJarDir = mfilename('fullpath');

folderSlash = '\';
if isunix
    folderSlash = '/';
end

juicerJarDirLevels = strfind(juicerJarDir,folderSlash);
juicerJarDir = [juicerJarDir(1:juicerJarDirLevels(end)),'juicer_tools.jar'];

% juicerJarDir = strrep(juicerJarDir,folderSlash,[folderSlash,folderSlash]);
% juicerJarDir = strrep(juicerJarDir,' ',[folderSlash,' ']);

% run juicer Dump
[status,cmdout] = system(sprintf('java -jar "%s" dump %s %s %s %s %s %s %s %s',...
    juicerJarDir,norm2d,norm1d,fn,chr1,chr2,bpFrag,binSize,fnOut));

% report an error if status~=0
if status~=0
    error(cmdout)
end

%format data
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


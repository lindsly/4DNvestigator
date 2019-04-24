function [normVec] = juicerToolsDumpNorm(norm1d,fn,chr,bpFrag,binSize,fnOut)
%juicerToolsDumpNorm MATLAB version of juicertools "dump norm" to read .hic
%   juicerToolsDumpNorm reads in .hic files and outputs the "norm1d"
%   normalization vector of chromosome "chr" at "binSize" resolution to the
%   file "fnOut".
%
%   Reference:
%   https://github.com/aidenlab/juicer/wiki/Data-Extraction#Examples
%
%   Input
%   norm1d:     1D normalization [NONE/VC/VC_SQRT/KR]
%   fn:         HicFile(s) location
%   chr1:       Chromosome # (eg 1-22,X,Y in human)
%   bpFrag:     Bin units [BP/FRAG] (FRAG dependant on RE)
%   binSize:    Bin size (ie 1E5 for 100kb resolution)
%   fnOut:      Temporary file name to output. Recommended that this not
%               input this, and let it be default (this will create a
%               temporary .txt file in your wd, which will be deleted
%               automatically)
%
%   Output
%   normVec: 	normalization vector for given Hi-C parameters 
%
%   Scott Ronquist, scotronq@umich.edu. 1/22/19

%% set default parameters
if ~exist('fnOut','var')||isempty(fnOut)
    fnOut = [tempdir, 'juicer_temp.txt'];
end

%% num2str if necessary
if isnumeric(chr);chr=num2str(chr);end
if isnumeric(binSize);binSize=num2str(binSize);end

%% format juicer tools input info
% find 4DNvestigator juicertools path
juicerJarDir = mfilename('fullpath');

juicerJarDirLevels = strfind(juicerJarDir,filesep);
juicerJarDir = [juicerJarDir(1:juicerJarDirLevels(end)),'juicer_tools.jar'];

%% run juicer Dump - Matrix
status = 1;
attempts = 1;
while status
    [status,cmdout] = system(sprintf('java -jar "%s" dump norm %s %s %s %s %s %s',...
        juicerJarDir,norm1d,fn,chr,bpFrag,binSize,fnOut));
    
    % report an error if status~=0, and we've tried twice
    if status~=0
        
        % try again if time out error
        if contains(cmdout,'Read timed out') && attempts < 3
            attempts = attempts+1;
        else
            error(cmdout)
        end
    end
end

%% format normalization vector
temp = readtable(fnOut,'Delimiter','\t');
normVec = temp{:,1};
if strcmp(fnOut,[tempdir, 'juicer_temp.txt'])
    delete(fnOut)
end

end


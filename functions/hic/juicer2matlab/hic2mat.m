function [H,hicHeader] = hic2mat(norm3d,norm1d,fn,chr1,chr2,bpFrag,binSize,intraFlag,headerFlag)
%hic2mat converts .hic files to MATLAB matrices
%
%   Reference:
%   https://github.com/theaidenlab/juicer/wiki/Data-Extraction
%
%   Inputs
%   norm_2d:    2D normalization [observed/oe]
%   norm_1d:    1D normalization [NONE/VC/VC_SQRT/KR]
%   fn:         hicFile(s) location
%   chr1:       chromosome # (eg 1-22,X,Y in human)
%   chr2:       chromosome # (eg 1-22,X,Y in human)
%   bpFrag:     bin units [BP/FRAG] (FRAG dependant on RE)
%   binSize:    bin size (ie 1E5 for 100kb resolution)
%   fnOut:      temporary file name to output. Recommended that this not
%               input this, and let it be default (this will create a
%               temporary .txt file in your wd, which will be deleted
%               automatically)
%
%   Output
%   H:          NxN double array for Hi-C matrix
%
%   Scott Ronquist, scotronq@umich.edu. 1/22/19

%% set default parameters
if ~exist('intraFlag','var')||isempty(intraFlag); intraFlag=0; end
if ~exist('headerFlag','var')||isempty(headerFlag); headerFlag=0; end

%% run juicertools and format
[juicerOut,hicHeader] = juicerToolsDump(norm3d,norm1d,fn,chr1,chr2,bpFrag,binSize,[],headerFlag);
H = juicerDump2mat(juicerOut,intraFlag);

end


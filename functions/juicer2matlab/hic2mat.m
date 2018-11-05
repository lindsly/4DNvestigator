function [H] = hic2mat(norm2d,norm1d,fn,chr1,chr2,bpFrag,binSize,intraFlag)
%hic2mat converts .hic files to MATLAB matrices
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
%   Output:
%   H: NxN double array for Hi-C matrix
%
%   Example:
%   
%
%   Scott Ronquist, 8/28/18

if ~exist('intraFlag','var')||isempty(intraFlag); intraFlag=0; end

H = juicerDump2mat(juicerToolsDump(norm2d,norm1d,fn,chr1,chr2,bpFrag,binSize),intraFlag);

end


function [HTrim,badLocs] = hicTrim(H,numDiag,numSparse)
%hic_trim trims Hi-C data to remove poor mapping regions
%   inputs
%   H: Hi-C matrix to be trimed (NxN double; default: N/A)
%   
%   numDiag: Number of bins off the diagonal to consider. All elements
%   further than numDiag from the diagonal will be considered (integer;
%   default: 1)
%
%   numSparse: specifies which bins to consider sparse. If 0<=numSparse<1,
%   this specifies the threshold proportion for keeping a row (ie rows with
%   <= numSparse bins full are removed). if numSparse>1, this is the number
%   bins that must be observed per row to be kept(type; default: )
%
%   outputs
%   HTrim: the trimmed Hi-C matrix (MxM double)
%   badLocs: location of removed locs, relative to the original Hi-C (type)
%
%   example
%   [B] = function(A)
%
%   Scott Ronquist, 7/29/18

%% set defaults if not specified
if nargin<2;numDiag=1;end
if nargin<3;numSparse=0;end

% turn nan to zero
H(isnan(H)) = 0;

%% loop through samples and check for bad_locs
badLocs = zeros(size(H,3),size(H,1));
for i = 1:size(H,3)
    % make sure symmetric
    HTemp = (H(:,:,i)+H(:,:,i)')./2;
    
    % remove diagonal
    HTemp = triu(HTemp-(triu(HTemp)-triu(HTemp,numDiag)));
    HTemp = HTemp+HTemp';
    
    % determine where matrix is zero
    badLocs(i,:) = sum(HTemp)==0;
    
    % remove if sparse
    if numSparse < 1
        temp = sum(logical(HTemp))./size(HTemp,1) < numSparse;
        badLocs(i,:) = logical(badLocs(i,:)+temp);
    else
        temp = sum(logical(HTemp)) < numSparse;
        badLocs(i,:) = logical(badLocs(i,:)+temp);
    end
    
end

%% trim matrix
%determine where there is a bad_loc in any sample and remove
badLocs = logical(sum(badLocs,1));

HTrim = H(~badLocs,~badLocs,:);

end


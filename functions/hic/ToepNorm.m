function [NX,zero_locs,Tp] = ToepNorm(X,zero_locs)
% function: nomalize the matrix by dividing the diagonal and each paralell
%   of the diagonal by its mean
% Inpust:
%     X - matrix to be normalized, can be either square or not
%     zero_locs - logical vector denoting regions to exclude from "expected
%     mean" calculation, unmappable regions
%
% Scott Ronquist, March 2018
% scotronq@umich.edu

if nargin < 2
    %     zero_locs = sum(X)==0;
    temp = X; temp = temp-diag(diag(temp));
    zero_locs = sum(temp)==0;
end
X_mappable = ones(size(X));
X_mappable(zero_locs,:) = 0;
X_mappable(:,zero_locs) = 0;

if isempty(X)
    NX = [];
    Tp = [];
    return;
end

% Get size informaiton
[m,n] = size(X);
ms = min(m,n);
mxs = max(m,n);

% Diagonal summation
ds = sumDiag(X);
% Number of elements
Ne = sumDiag(X_mappable>0); 
Ne = Ne(end:-1:1);
% Diagonal mean value
mds = ds(end:-1:1)./Ne;

% Normalization matrix
Tp = toeplitz(mds(m:-1:1),mds(m:end));
% Nomralization
NX = X./Tp;

NX(isinf(NX))=0;
NX(isnan(NX))=0;
end

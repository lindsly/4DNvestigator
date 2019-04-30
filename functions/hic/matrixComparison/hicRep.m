function [ps] = hicRep(X,binSize)
%hicRep performs HiCRep matrix comparison
%   https://genome.cshlp.org/content/early/2017/08/30/gr.220640.117
%   
%   Inputs
%   X:          NxNx2 raw Hi-C matrix	
%   binSize:	Hi-C bin size to determine smoothing parameter (bp)
%   
%   Outputs
%   ps:         stratum adjusted correlation coefficient (SCC)
%   
%   Version 1.0 (04/29/19)
%   Written by: Scott Ronquist
%   Contact: 	scotronq@umich.edu
%   Created: 	04/29/19
%   
%   Revision History:
%   v1.0 (04/29/19)
%   * hicRep.m created

%% set up
X(isnan(X)) = 0;
N = length(X);
T = size(X,3); % should always be 2 (pairwise comparison)

% select h (smoothing parameter size, based on input Hi-c bin size)
hMat = [20  10E3;...
        11  25E3;...
        5   40E3;...
        3   100E3;...
        1   500E3;...
        0   1E6];
h = hMat(hMat(:,2)==binSize,1);

% smooth the matrix
hSmoothMat = ones(2*h+1)/(2*h+1)^2;
XSmooth = imfilter(X,hSmoothMat,'replicate');

% remove counts > 5Mb from diagonal
if N*binSize > 5E6
    remMat = triu(ones(N),round(5E6*(1/binSize)));
    remMat = remMat+remMat';
    XSmooth(find(repmat(remMat,[1 1 T]))) = NaN;
end

% Stratum-adjusted correlation coefficient (SCC)
kTotal = sum(~isnan(XSmooth(1,:,1)));
Nk = zeros(kTotal,1);
r1k = zeros(kTotal,1);
r2k = zeros(kTotal,1);
pk = zeros(kTotal,1);

for iK = 1:kTotal
    
    % get stratum contacts
    Xk = diag(XSmooth(:,:,1),iK-1);
    Yk = diag(XSmooth(:,:,2),iK-1);
    Nk(iK) = length(Xk);
    
    % calculate stratum mean and variance
    r1k(iK) = mean(Xk.*Yk)-(mean(Xk)*mean(Yk));
    r2k(iK) = sqrt(var(Xk)*var(Yk));
    r2kNorm(iK) = sqrt(var(normalize(Xk,'range'))*var(normalize(Yk,'range')));
    
    % stratum-specific correlation ?k
    pk(iK) = r1k(iK)/r2k(iK);
end

% stratum-adjusted correlation coefficient (SCC)
ps = sum(Nk.*r2kNorm.*pk)/sum(Nk.*r2kNorm);

end

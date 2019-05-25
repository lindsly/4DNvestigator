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
n = size(X,1);
T = size(X,3); % should always be 2 (pairwise comparison)

% select h (smoothing parameter size, based on input Hi-c bin size)
hMat = [20  10E3;...
        11  25E3;...
        5   40E3;...
        5   50E3;...
        3   100E3;...
        1   500E3;...
        0   1E6];
h = hMat(hMat(:,2)==binSize,1);

% smooth the matrix
hSmoothMat = ones(2*h+1)/(2*h+1)^2;
XSmooth = imfilter(X,hSmoothMat);%,'replicate');

% manual smooth
XSmooth = zeros(size(X));
for iSamp = 1:T
    for i = 1:n
        for j = i:n
            XSmooth(i,j,iSamp) = sum(sum(X(max([1 i-h]):min([i+h n]),...
                max([1 j-h]):min([j+h n]),iSamp)))/(1+2*h)^2;
        end
    end
    XSmooth(:,:,iSamp) = XSmooth(:,:,iSamp)+triu(XSmooth(:,:,iSamp),1)';
end

% remove counts > 5Mb from diagonal
if n*binSize > 5E6
    remMat = triu(ones(n),round(5E6*(1/binSize)));
    remMat = remMat+remMat';
    XSmooth(find(repmat(remMat,[1 1 T]))) = NaN;
end

% Stratum-adjusted correlation coefficient (SCC)
kTotal = sum(~isnan(XSmooth(1,:,1)));
Nk = zeros(kTotal,1);
r1k = zeros(kTotal,1);
r2k = zeros(kTotal,1);
r2kNorm = zeros(kTotal,1);
pk = zeros(kTotal,1);

for iK = 1:kTotal
    
    % get stratum contacts
    Xk = diag(XSmooth(:,:,1),iK-1);
    Yk = diag(XSmooth(:,:,2),iK-1);
    Nk(iK) = length(Xk);
    
    % calculate stratum mean and variance
    r1k(iK) = mean(Xk.*Yk)-(mean(Xk)*mean(Yk));
    r2k(iK) = sqrt(var(Xk)*var(Yk));
    % r2kNorm(iK) = sqrt(var(normalize(Xk,'range'))*var(normalize(Yk,'range')));
    
    % https://github.com/qunhualilab/hicrep/blob/master/R/vstran.R
    dataSorted = sort(Xk); [~, rnkXk] = ismember(Xk,dataSorted);
    dataSorted = sort(Yk); [~, rnkYk] = ismember(Yk,dataSorted);
    r2kNorm(iK) = sqrt(var(rnkXk/Nk(iK))*var(rnkYk/Nk(iK)));
    
    % stratum-specific correlation ?k
    pk(iK) = r1k(iK)/r2k(iK);
    
    %https://github.com/qunhualilab/hicrep/blob/master/R/get.scc.R
    pk(iK) = corr(Xk,Yk);
end

% stratum-adjusted correlation coefficient (SCC)
ps = nansum(Nk.*r2kNorm.*pk)/nansum(Nk.*r2kNorm);

end

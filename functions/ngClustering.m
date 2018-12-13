function [idx,clusterMap] = ngClustering(dataIn,K,sigma,plotFlag)
%ngClustering performs the Andrew Ng clustering algorithm
%   https://papers.nips.cc/paper/2092-on-spectral-clustering-analysis-and-an-algorithm.pdf
%
%   Input:
%   dataIn: set of points "S" or affinity matrix "A". Function can take
%   either input, determines this automatically
%   K: number of clusters for kmeans
%   sigma: scaling parameter
%   plotflag: logical for whether to plot the clustered results
%
%   Output:
%   idx: clustering index results
%   clusterMap: clusterin results overlayed to affinity matrix
%
%   Scott Ronquist, scotronq@umich.edu. 12/11/18

%% set default parameterss
constructAFlag = 1;
if issymmetric(dataIn);constructAFlag = 0;end
if ~exist('K','var')||isempty(K);K=5;end
if ~exist('sigma','var')||isempty(sigma);sigma=1;end
if ~exist('plotFlag','var')||isempty(plotFlag);plotFlag=0;end

%% Ng Algorithm
% 1. Form the affinity matrix A, and Aii = 0
if constructAFlag
    A = tril(ones(size(s,1)),-1);
    A(logical(A))= pdist(dataIn);
    A = A + A.';
    A = exp((-A.^2)/sigma);
else
    A = dataIn;
end
A = A - diag(diag(A));

% 2. Define D to be the diagonal matrix and construct the matrix
% L = D^(-l/2)*A*D^(-l/2)
D = diag(sum(A));
D_inv = D^(-1/2);
L = D_inv*A*D_inv;

% 3. Find xl ,x2, ..., xk, the k largest eigenvectors of L, and form the
% matrix X = [xl ,x2, ..., xk];
[eigVec,~] = eigs(L,K,'largestreal');%'largestreal');largestabs
X = eigVec;

% 4. Form the matrix Y from X by renormalizing each of X's rows to have
% unit lengths
Y = diag(1./sqrt(sum(X.*X,2)))*X;

% 5. Treating each row of Y as a point in Rk , cluster them into k clusters
% via K-means
[idx, ~] = kmeans(Y,K,'start',eye(K));

%% create clusterMap
clusterMap = zeros(size(X));
for i = 1 : K
    tempIdx = find(idx ==i);
    clusterMap(tempIdx,tempIdx)=i;
end
clusterMap = clusterMap*max(A(:));

%% plot
if plotFlag
    figure
    subplot(1,2,1)
    imagesc(A), axis square
    title('input matrix')
    
    subplot(1,2,2)
    imagesc(A), axis square, hold on
    h = imagesc(clusterMap);
    set( h, 'AlphaData', .7 );
    title('clustered matrix')
    
    set(get(gcf,'children'),'fontsize',15,'linewidth',2)
    linkaxes(get(gcf,'children'))
end

end


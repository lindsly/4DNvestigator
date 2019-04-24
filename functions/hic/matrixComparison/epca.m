function [epcaScore,epcaMat] = epca(A)
%EPCA performs E-PCA
%   
%   Input
%   A: N x N x M contact matrix
%
%   Output
%   epcaScore:  PCA scores
%   epcaMat:    PCA scores mapped to N x N
%
%   Reference:
%   - Characterizing the 3D structure and dynamics of chromosomes and
%   proteins in a common contact matrix framework.
%   https://academic.oup.com/nar/article/46/16/8143/5051111
%
%   - website for "CAMERRA" (used for IPCA and EPCA)
%   http://shenlab.utk.edu/camerra.html
%
%
%   Version 1.0 (4/14/19)
%   Written by: Scott Ronquist
%   Contact: scotronq@umich.edu
%   Contributors:
%   Created: 4/14/19
%   Revision History:

%% Set-up
% create random data for testing
if 1==0
    N = 10;
    T = 4;
    A = zeros(N, N, T);
    for i = 1:T
        temp = rand(N);
        A(:,:,i) = triu(temp)+triu(temp,1)';
    end
end

% get input data parameters
N = size(A,1);
T = size(A,3);

% vectorize input matrices
Avec = zeros(N^2,T);
for i = 1:T
    temp = A(:,:,i);
    Avec(:,i) = temp(:);
end

% calculate covariance matrix and pca
% *NOTE: can be optimized in the future for symmetric matrices, only need
% (N choose 2)+N, not N^2
covA = cov(Avec');
[~, epcaScore] = pca(covA);

% map points back to N x N space
epcaMat = zeros(N,N,T);
for i = 1:T
    epcaMat(:,:,i) = reshape(epcaScore(:,i),N,N);
end

end


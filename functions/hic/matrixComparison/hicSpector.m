function [Sd,Q] = hicSpector(A,B,r)
%hicSpector method for Hi-C comparison
%   hicSpector computes the spectral decomposition of the normalized
%   laplacian of two different Hi-C matrices, and provides a score from
%   [0,1] for "reproducibility", ie how similar the matrices are. 0 means
%   matrices are highly similar, 1 means matrices are highly different.
%
%   recommended input is raw contact matrices, 40kb bin, intra-chr
%
%   A: Hi-C matrix from sample A
%   B: Hi-C matrix from sample B
%   r: number of eigenvectors to compare (default:20)
%
%   Sd: eigenvector difference summation
%   Q: reproducibilty score [0 1]
%
%   Reference:
%   Yan, Koon-Kiu, et al. "HiC-spector: a matrix library for spectral and
%   reproducibility analysis of Hi-C contact maps." Bioinformatics 33.14
%   (2017): 2199-2201.
%
%   Version 1.0 (4/14/19)
%   Written by: Scott Ronquist
%   Contact: scotronq@umich.edu
%   Contributors:
%   Created: 4/14/19
%   Revision History:

%% set-up
% default number of eigenvectors
if nargin<3; r=min([20, size(A,1)]); end

% if matrix size >= r, r = matrix size-1
if r >= size(A,1)
    r = size(A,1)-1;
end

% calculate the laplacian from Hi-C
LnA = mat2lap(A);
LnB = mat2lap(B);

% get the top r eigenvectors
[LnAeigvec,~] = eigs(LnA,r+1,'SM');
[LnBeigvec,~] = eigs(LnB,r+1,'SM');

% distance metric
Sd = 0;
for i = 1:r
    % eigenvector direction is arbitrary to take min(vA-vB,vA+vB)
    tempDist = min([norm(LnAeigvec(:,i+1)-LnBeigvec(:,i+1)),...
        norm(LnAeigvec(:,i+1)+LnBeigvec(:,i+1))]);
    Sd = Sd+tempDist;
end

% reproducibility score (bounded [0 1])
Q=(1-(1/r)*(Sd/sqrt(2)));

    % function to calculate the Laplcian matrix
    function [L] = mat2lap(A)
        D = diag(sum(A))^(-.5);
        L = (D^(-.5))*(D-A)*(D^(-.5));
    end
end



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
%   Scott Ronquist, 6/27/18

if nargin<3;r=20;end

Ln_A = hic_laplacian_fdv(A);
Ln_B = hic_laplacian_fdv(B);

[Ln_A_v,~] = eigs(Ln_A,r+1,'SM');
[Ln_B_v,~] = eigs(Ln_B,r+1,'SM');

Sd = 0;
for i = 1:r
    % eigenvector direction is arbitrary to take min(vA-vB,vA+vB)
    tempDist = min([norm(Ln_A_v(:,i+1)-Ln_B_v(:,i+1)),...
        norm(Ln_A_v(:,i+1)+Ln_B_v(:,i+1))]);
    Sd = Sd+tempDist;
end

Q=(1-(1/r)*(Sd/sqrt(2)));
end


function [Ln,Fdv,L,FdNum] = hicLaplacianFdv(A)
%hicLaplacianFdv calculates the graph laplacian from the matrix A
%   inputs
%   A: Hi-C contact matrix, typically normalized (NxN double; default: N/A)
%
%   outputs
%   Ln: laplacian matrix, normalized (NxN double)
%   Fdv: Fielder vector (Nx1 double)
%   L: laplacian matrix, unnormalized (NxN double)
%
%   example
%   [B] = hicLaplacianFdv(A)
%
%   Scott Ronquist, 7/29/18

%% calculate normalized laplacian
D = diag(sum(A));
L = D-A;
Dn = diag(diag(D).^(-1/2));
Ln = Dn*(D-A)*Dn;

%% calculate eigenvalue and eigenvector
try
    [F,D] = eigs(Ln,2,'SM');
    Fdv = F(:,2);
    FdNum = D(2,2);
catch
    fprintf('eigs ERROR: First input matrix is singular. trying full eig \n')
    [V,D] = eig(Ln);
    [~,I] = sort(diag(D),'ascend');
    Fdv = V(:,I(2));
    FdNum = D(I(2),I(2));
end
end
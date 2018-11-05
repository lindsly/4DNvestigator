function [vnEntropy,AInputEigs] = hicVnEntropy(AInput,numEigs,normEigs,preProcess)
%hicVnEntropy computes the Von Neumann graph entropy of the matrix A
%
%   inputs
%   AInput: Input matrix for Von Neumann entropy calculation (NxN double; default: N/A)
%
%   numEigs: number of eigenvalues to consider for compuation (integer;
%   default: min(size(AInput,1),20))
%
%   normEigs: normalize eigenvalues flag ([0,1]; default: 1)
%
%   preProcess: arguments for preprocessing methods (string; default: 'none')
%
%   outputs
%   vnEntropy: Von Neumann Entropy of the input matrix (double)
%
%   AInputEigs: eigenvalues of the input matrix (numEigsx1 double)
%
%
%   example
%   [B] = function(A)
%
%   Scott Ronquist, 7/29/18

if nargin < 4;preProcess='none';end
if nargin < 3;normEigs=1;end
if nargin < 2;numEigs=min(size(AInput,1),20);end

for i = 1:size(AInput,3)
    %% preprocess
    switch preProcess
        case 'none'
            A = AInput(:,:,i);
        case 'laplacian'
            A = hic_laplacian_fdv(AInput(:,:,i));
        case 'corr'
            A = corr(AInput(:,:,i));
            A(isnan(A))=0;
    end
    
    %% eigen decomposition
    [V,D] = eigs(A,numEigs);
    AInputEigs(:,i) = diag(D);
    
    %% normalize eigenvalues
    if normEigs == 1
        AInputEigs(:,i) = AInputEigs(:,i)./sum(AInputEigs(:,i));
    end
    
    %% Von Neumann Entropy
    vnEntropy(i) = -sum(AInputEigs(:,i).*log(AInputEigs(:,i)));
end

end


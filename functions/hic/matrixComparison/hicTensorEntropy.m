function [entropyMatrix] = hicTensorEntropy(AInput,thresh)
%hicVnEntropy computes the tensor entropy of the input time series Hi-C matrices
%
%   Inputs
%   AInput:     Input matrices for entropy calculation 
%               (NxNxM double; default: N/A)
%
%   Outputs
%   entropy: Tensor entropy of the input matrices (double)

%% Set Parameters
    k = 3;
    if ~exist('thresh','var')||isempty(thresh);thresh = .95;end

%% Calculate multicorrelation and tensor entropy
    entropyMatrix = zeros(1, size(AInput, 3));

    for i = 1:size(AInput, 3)
        % Compute multicorrelations using input matrix
        corrMatrix = corr(log2(AInput(:, :, i)));
        n = size(corrMatrix, 1);
        edges = nchoosek(1:n, k);
        multiCorr = zeros(size(edges, 1), 1);
        for j = 1:size(edges, 1)
            R =corrMatrix(edges(j, :), edges(j, :));
            multiCorr(j) = (1-det(R))^(1/2);
        end
        createdEdges = edges(multiCorr>=thresh, :);

        % Construct the adjacency and degree tensors
        adjTensor = zeros(n*ones(1, k));
        degTensor = zeros(n*ones(1, k));
        for j = 1:size(createdEdges, 1)
            allPerms = perms(createdEdges(j, :));
            for jj = 1:size(allPerms, 1)
                % When adjusting k, add terms here accordingly
                adjTensor(allPerms(jj, 1), allPerms(jj, 2), allPerms(jj, 3)) = 1/factorial(k-1);
            end
        end
        for j = 1:n
            % When adjusting k, add terms here accordingly
            degTensor(j, j, j) = sum(adjTensor(j, :), 'all');
        end
        
        % Compute the tensor entropy
        lapTensor = degTensor-adjTensor;
        reLapMatrix = reshape(lapTensor, n^(k-1), n);
        genSingValues = svd(reLapMatrix, 'econ');
        genSingValues = genSingValues/sum(genSingValues);
        genSingValues = genSingValues(genSingValues>0);
        entropyMatrix(1, i) = -sum(genSingValues.*log(genSingValues));
    end

end
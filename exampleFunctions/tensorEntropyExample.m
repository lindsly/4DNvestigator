tic
load('myodEntropy.mat');

k = 3;

combData = cat(3, chr14_fib, chr14_myod);
eps = 0.95;

entropyMatrix = zeros(2, size(combData, 3));
for i = 1:size(combData, 3)
    % Network Entropy
    adjMatrix = combData(:, :, i);
    adjMatrix = adjMatrix - diag(diag(adjMatrix));
    adjMatrix = adjMatrix>=0.95;
    lapMatrix = diag(sum(adjMatrix, 1))-adjMatrix;
    eigValues = eig(lapMatrix);
    eigValues = eigValues/sum(eigValues);
    eigValues = eigValues(eigValues>0);
    entropyMatrix(1, i) = -sum(eigValues.*log(eigValues));
    
    
    % Tensor Entropy 
    corrMatrix = corr(log2(combData(:, :, i)));
    n = size(corrMatrix, 1);
    edges = nchoosek(1:n, k);
    multiCorr = zeros(size(edges, 1), 1);
    for j = 1:size(edges, 1)
        R =corrMatrix(edges(j, :), edges(j, :));
        multiCorr(j) = (1-det(R))^(1/2);
    end
    createdEdges = edges(multiCorr>=eps, :);
    
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
    lapTensor = degTensor-adjTensor;
    reLapMatrix = reshape(lapTensor, n, n^(k-1));
    genSingValues = svd(reLapMatrix, 'econ');
    genSingValues = genSingValues/sum(genSingValues);
    genSingValues = genSingValues(genSingValues>0);
    entropyMatrix(2, i) = -sum(genSingValues.*log(genSingValues));
end

figure
plot(1:8, entropyMatrix(1,9:16)) % Reprogramming
hold on
plot(1:8,entropyMatrix(1,1:8)) % Proliferation
set(gca,'YLim',[3 5])

figure
plot(1:8, entropyMatrix(2,9:16)) % Reprogramming
hold on
plot(1:8,entropyMatrix(2,1:8)) % Proliferation
set(gca,'YLim',[3 5])

toc

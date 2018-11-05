% Preprocessing of the adjacency matrix A of a graph to produce
% a symmetric matrix, VSsims, of the cosine similarities of
% V, the right singular vector of A, and S, the singular values
% of A. VSsims is then thresholded into a sparse symmetric
% matrix, VSsims t, by eliminating those values of VSsims which
% are effectively 0. This implementation uses a threshold of 0.2
function [VS, VSsimsT] = preprocess2(A)
t = 0.2;
k = 30;
[numTerms,numDocuments] = size(A);
[U,S,V] = svds(A,k);
VS = V*S;
for j = 1:numDocuments
    VS(j,:) = VS(j,:) / norm(VS(j,:),2);
end
VSsims = VS*VS';
VSsimsThresh = VSsims > .2;
VSsimsNew = VSsims.*VSsimsThresh;
VSsimsT = sparse(VSsimsNew);
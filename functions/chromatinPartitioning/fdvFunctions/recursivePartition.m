% Overarching method face to allow for various algorithms to be
% used in partitioning. The desired algorithm is specified via
% the 3rd paramter, a String denoting the name of the algorithm
% to use. This method then calls the appropriate partitioning
% method where all necessary variables are assigned.
function [partitions, index] = recursivePartition(A,n,criteria)
[~, A] = preprocess2(A);
if strcmpi(criteria, 'normalized')
    [partitions, index] = recursiveNormalized(A,n);
elseif strcmpi(criteria, 'ratio')
    [partitions, index] = recursiveRatio(A,n);
elseif strcmpi(criteria, 'min')
    [partitions, index] = recursiveMinCut(A,n);
elseif strcmpi(criteria, 'minmax')
    [partitions, index] = recursiveMinMax(A,n);
elseif strcmpi(criteria, 'bcut')
    [partitions, index] = recursiveBCut(A,n);
else
    error('Unknown partitioning criteria');
end


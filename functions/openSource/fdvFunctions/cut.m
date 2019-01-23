% Function to calculate cut(A,B), ie the set of edges between two
% disjoint sets A,B both in the graph I, where A+B=I, the
% complete graph.
% @precondition: subset A is the first portion of the graph,
% the graph, sortedI, and subset B is the distinct second
% portion of sortedI.
% A and B are square symmetric adjacency matrices.
function edgeTotal = cut(sortedI,sizeA,sizeI)
discard = sortedI(sizeA+1:sizeI, 1:sizeA);
edgeTotal = sum(sum(discard));


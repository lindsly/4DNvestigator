function [effRes] = effectiveResistance(A)
%effectiveResistance calculates the effective resistance of an adjacency
%   matrix
%
%   Input
%   A: adjacency matrix
%
%   Output
%   effRes = effective resistance
%
%   Example
%   s = [1 1 1 2 2 3];
%   t = [2 3 4 3 4 4];
%   G = graph(s,t);
%   A = full(adjacency(G));
%   effectiveResistance(A)
%
%   Resources
%   https://beta.vu.nl/nl/Images/werkstuk-meer_tcm235-280356.pdf
%   https://www.nas.ewi.tudelft.nl/people/Piet/papers/LAA_2011_EffectiveResistance.pdf
%
%   Scott Ronquist, scotronq@umich.edu. 12/11/18

%% calculate effective resistance
D = diag(sum(A));
L = D-A;
Ln = (D^(-1/2))*L*(D^(-1/2));       % normalized laplacian, not sure if we should use this yet...
lambda = eig(L);
lambda = sort(lambda,'ascend');
effRes = length(A)*sum(1./lambda(2:end));

end


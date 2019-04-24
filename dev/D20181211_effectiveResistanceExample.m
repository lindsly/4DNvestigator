% Goal: effective resistance examples from 
% https://beta.vu.nl/nl/Images/werkstuk-meer_tcm235-280356.pdf
% https://www.nas.ewi.tudelft.nl/people/Piet/papers/LAA_2011_EffectiveResistance.pdf

%% a
clear
close all
s = [1 1 1 2 2 3];
t = [2 3 4 3 4 4];
G = graph(s,t);
A = full(adjacency(G));
effectiveResistance(A)

%% b
s = [1 1 1];
t = [2 3 4];
G = graph(s,t);
A = full(adjacency(G));
effectiveResistance(A)


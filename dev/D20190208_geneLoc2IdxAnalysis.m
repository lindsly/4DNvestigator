% goal: script to convert gene location information to corresponding
% trimmed index

%% set up
clear
close all

%% create sample data
genePos = [ 100 900;...
            85 702];
trimLocs = zeros(100,1);
trimLocs(2) = 1;
trimLocs(70:71) = 1;

%% run geneLoc2Idx
[genePosNew] = geneLoc2Idx(genePos,10,trimLocs);
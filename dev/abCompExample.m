% AB compartmentalization example with different fiedler cut methods
% https://www.sciencedirect.com/science/article/pii/S0024379506003454

clear
close all

%% create example data
rng(1);
numNodes = 100;

%% create adjacency matrix
a = rand(numNodes);
if 1==1
    a = a.*(rand(numNodes)>.7);
end
a = triu(a)+triu(a,1)';
figure, plot(graph(a),'EdgeAlpha',.2,'Layout','force')
figure, imagesc(a)

% load barbellgraph.mat
% a = full(A);

%% 

%% run through AB comp analysis
fdvSplitMethods = {'gap'};%{'bisection','ratio','sign','gap'};
for i = 1:length(fdvSplitMethods)
    [ABcomp,abGroup] = hicABcomp(a,[],[],[],[],1,fdvSplitMethods{i});
    
    % graph plot
    figure
    p = plot(graph(a),'Marker','.','EdgeAlpha',.2,'Layout','force','MarkerSize',10);
    highlight(p,find(abGroup==1),'NodeColor','r') % subgraph A
    highlight(p,find(abGroup==2),'NodeColor','k') % subgraph B
end

%% EXTRA
[ABcomp,abGroup] = hicABcomp(a,[],[],[],[],1,'ratio');



% Goal: show Centrality + dimension reduction for real data
%   Hi-C + RNA-seq, MYOD; 2 time points, 1 chr
%
% Hi-C data
% MYOD 2 time points Hi-C RNa-seq
%
% feature extraction (centrality)
% same space (normalize features)
% dimension reduction (Lapcian eigenmaps)
%
%   Scott Ronquist, scotronq@umich.edu. 3/28/19

%% Set up
clear
close all

addpath('data')
addpath('functions')
addpath(genpath('drtoolbox'))

%% Parameters
graphWeighted=1;
binInfo=[];
norm=0;

%% load data
load('MyoD_kb_HiC_rpm.mat', 'C14')
load('Sijia_100kb_bin_RNAseq.mat', 'RNAseq_bin_myoD','list_gene_name_bin_myoD')

% format data
H = C14(:,:,1:2);
R = RNAseq_bin_myoD{14}(:,1:2);
binNames=list_gene_name_bin_myoD{14};

% normalize Hi-C
for iTp = 1:size(H,3)
    H(:,:,iTp) = ToepNorm(H(:,:,iTp));% figure, imagesc(log(H(:,:,1)))
end

% trim Hi-C (remove disconnected nodes)
thresh = 0;
badLocs = any(squeeze(sum(H,1)) <= thresh, 2);
H(badLocs,:,:) = [];
H(:,badLocs,:) = [];
R(badLocs,:) = [];
binNames(badLocs) = [];

% visualize raw data
figure
subplot(6,1,1)
bar(log(R(:,1)+1))

subplot(6,1,2:6)
imagesc(log(H(:,:,1)))

%% extract features (centrality)
for tp  = 1:size(H,3)
    G = graph(H(:,:,tp),'OmitSelfLoops');
    
    % from Sijia, inverse weights for closeness and betweeness
    Ht_root = nthroot(H(:,:,tp),2);
    Ct = 1./Ht_root;
    Ct(isnan(Ct)) = 0;  Ct(isinf(Ct)) = 0;
    G_inv = graph(Ct,'OmitSelfLoops');
    
    % Extract features
    if graphWeighted % centrality measure cost/importance is edge weighted
        features(:,:,tp) = [centrality(G,'degree','Importance',G.Edges.Weight),...
            centrality(G,'closeness','Cost',G_inv.Edges.Weight),...
            centrality(G,'betweenness','Cost',G_inv.Edges.Weight),...
            centrality(G,'eigenvector','Importance',G.Edges.Weight)];
    else
        features(:,:,tp) = [centrality(G,'degree'),centrality(G,'closeness'),...
            centrality(G,'betweenness'),centrality(G,'eigenvector')];
    end
end

% add RNA-seq and normalize features (sam length vector)
features = cat(2,features,reshape(R,size(R,1),1,size(R,2)));
Xnorm = FeatureNorm2(features);

% stack time points
XnormStacked = [];
labels = [];
for i = 1:size(Xnorm,3)
    XnormStacked = [XnormStacked;Xnorm(:,:,i)];
    labels = [labels;ones(size(Xnorm(:,:,i),1),1)*i];
end
uniqLabels = unique(labels);

%% dimension reduction (Lapcian eigenmaps)
dimRedType = 'lapEigen';

switch dimRedType
    case 'pca'
        [mappedX, mapping] = pca(XnormStacked);
    case 'tsne'
        mappedX = tsne(XnormStacked);
    case 'lapEigen'
        [mappedX, mapping] = laplacian_eigen(XnormStacked);
    otherwise
        error('Bioinf547DimReduc: invalid dimension reduction type "%s"',dimRedType)
end

%% figure
% get even spaced colors
colormapJet = jet;
colors = colormapJet(ceil(linspace(1,size(colormapJet,1),length(uniqLabels))),:);

% visualize data
figure, hold on
for iLabel = 1:length(uniqLabels)
    labelLoc = labels==uniqLabels(iLabel);
    plot(mappedX(labelLoc,1),mappedX(labelLoc,2),...
        '.','color',colors(iLabel,:),'markersize',20)
end
legend({'before MYOD1','after MYOD1'})

% format figure
box on
title(sprintf('Data: "%s", Dim reduc: "%s"','Hi-C',dimRedType), 'Interpreter', 'none')
set(get(gcf,'children'),'fontsize',15,'linewidth',2)

% add line between points
plot([mappedX(1:883,1) mappedX(884:end,1)]',...
    [mappedX(1:883,2) mappedX(884:end,2)]','-')

% add gene names to top N shift points
XChange = zeros(883,1);
for i = 1:883
    XChange(i) = sqrt((mappedX(i,1)-mappedX(i+883,1))^2 +...
        (mappedX(i,2)-mappedX(i+883,2))^2);
end

[B,idx] = sort(XChange,'descend');

% add gene name text
for i = 1:10
    if isempty(binNames{idx(i)});continue;end
    text(mappedX(idx(i),1),mappedX(idx(i),2),binNames{idx(i)})
end


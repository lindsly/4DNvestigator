function [features,score,genes_change] = sfAnalysis_temp(hic,feature1D,binNames,binInfo,graphWeighted,dimReduc,topEllipseFrac,sfType,featureType)
%sfAnalysis analyzes Hi-C and RNA-seq through centrality and PCA
%   sfAnalysis extracts centrality features from Hi-C, concatenates these
%   features with RNA-seq, then determines a low dimensional projection to
%   extract regions that change significantly
%
%   Input
%   hic:            Hi-C matrix, typically normalized (NxNxM double;
%                   default: N/A)
%   rnaSeq:         Rna-seq values, (NxM double; default: N/A)
%   binNames:       Names associated with Hi-C and RNA-seq bins (string or
%                   cell array; default: empty cell)
%   norm:           Normalize Hi-C flag, observed over expected ([0,1]; default: 0)
%   binInfo:        Bin chr information necessary for normalization (array;
%                   default: [])
%   graphWeighted:  Use graph weighted centrality analysis ([0,1]; default: 1)
%   dimReduc:       Type of dimension reduction (PCA, LapEigen)
%   topEllipseFrac: The top fraction of bins that are fit with an ellipse
%
%   Output
%   features:       Concatenated feature array (Nx3 double)
%   score:          Location of bin in low dimensional projection (Nx3 double)
%
%   Version 1.1 (4/26/19)
%   Written by: Scott Ronquist
%   Contact:    scotronq@umich.edu
%   Created:    1/23/19
%   
%   Revision History:
%   v1.0 (1/23/19)
%   * sfAnalysis.m created
%   v1.1 (4/26/19)
%   * code commenting
%   * fixed ellipse fitting
%   * added dimReduc option

%% Set default parameters
if ~exist('dimReduc','var')||isempty(dimReduc);dimReduc='pca';end
if ~exist('graphWeighted','var')||isempty(graphWeighted);graphWeighted=1;end
if ~exist('binInfo','var')||isempty(binInfo);binInfo=[];end
% if ~exist('norm','var')||isempty(norm);norm=0;end
if ~exist('binNames','var')||isempty(binNames);binNames=cellstr(num2str([1:length(feature1D)]'));end
if ~exist('topEllipseFrac','var')||isempty(topEllipseFrac);topEllipseFrac=.1;end

load('ncbi_gene_info.mat')
color_opts = color_options;

%% Normalize, remove centromere and extract centrality (depreciated)
% if norm
%     error(['norm = 1 is currently a deprecated option, normalize Hi-C',...
%         'matrices prior to input'])
%     hic = norm_hic_bins(hic,binInfo);
% end

%% Extract Features
% extract centrality
for tp  = 1:size(hic,3)
    G = graph(hic(:,:,tp),'OmitSelfLoops');
    
    % Inverse weights for closeness and betweeness (MATLAB weight definition)
    HtRoot = nthroot(hic(:,:,tp),2);
    Ct = 1./HtRoot;
    Ct(isnan(Ct)) = 0;  Ct(isinf(Ct)) = 0;
    GInv = graph(Ct,'OmitSelfLoops');
    
    % Extract features
    if graphWeighted % centrality measure cost/importance is edge weighted
        features(:,:,tp) = [centrality(G,'degree','Importance',double(G.Edges.Weight)),...
            centrality(G,'closeness','Cost',double(GInv.Edges.Weight)),...
            centrality(G,'betweenness','Cost',double(GInv.Edges.Weight)),...
            centrality(G,'eigenvector','Importance',double(G.Edges.Weight))];
    else
        features(:,:,tp) = [centrality(G,'degree'),centrality(G,'closeness'),...
            centrality(G,'betweenness'),centrality(G,'eigenvector')];
    end
end

% use log scale to make better visualizations of RNA-seq
% adding .5 and 1 prevents negative RNA-seq values after log2 
if strcmp(featureType,'rna')
    feature1D = log2(feature1D+.5)+1;
end
% add extra feature to centrality measures (default is RNA-seq)
features = cat(2,features,reshape(feature1D,size(feature1D,1),1,size(feature1D,2)));

% normalize features
Xnorm = FeatureNorm2(features);

% stack time points to project to same low dimensional space
XnormStacked = [];
for i = 1:size(Xnorm,3)
    XnormStacked = [XnormStacked;Xnorm(:,:,i)];
end

if strcmp(sfType,'sfmatrix')
    %% Dimension reduction
    dimReduc = lower(dimReduc);
    dimReducLabel = [];
    switch dimReduc
        case 'pca'
            [~,score,~] = pca(XnormStacked);
            dimReducLabel = 'PCA';
        case 'lapeigen'
            knn = 30*size(hic,3);%30*size(hic,3);
            sigma = 100;
            score = laplacian_eigen(XnormStacked, 3, knn, sigma);
            dimReducLabel = 'Laplacian Eigenmaps';
        case 'tsne'
            score = tsne(XnormStacked,[],3);
            dimReducLabel = 't-SNE';
        case 'umap'
            [score,~,clust] = run_umap(XnormStacked);
            set(gcf,'Name','UMAP default plot')
            dimReducLabel = 'UMAP';
        otherwise
            error('please select a valid dimension reduction method: pca, lapeigen, tsne, or umap')
    end

    %% Figure output
    figure('Name',['Structure_Function_Feature_Space_', dimReducLabel],'position',[50 50 750 550]), hold on
    numSamples = size(Xnorm,3);
    numBins = size(Xnorm,1);
    scoreReshape = zeros(numBins,2,numSamples);
    for iS = 1:numSamples
        scoreReshape(:,:,iS) = score(numBins*(iS-1)+1:numBins*iS,1:2);
    end

    % get dist btw pts
    allDist = zeros(numBins,numSamples-1);
    for iS = 1:numSamples-1
        allDist(:,iS) = diag(pdist2(scoreReshape(:,1:2,iS), scoreReshape(:,1:2,iS+1)));
    end

    numPts_plot = size(scoreReshape,1);
    numPts_highlight = min(size(scoreReshape,1), 10);
    [sorted_allDist,labelLocs] = sort(sum(allDist,2),'descend');
%     sorted_allDist

    % plot features in low dimensional space
    colorScale = jet(numSamples);
    for iS = 1:numSamples
         scatter(scoreReshape(labelLocs(1:numPts_plot),1,iS), ...
             scoreReshape(labelLocs(1:numPts_plot),2,iS), ...
             15,color_opts(iS,:),'filled');
    end
    
    %%
    genes_change = binNames(labelLocs(1:numPts_highlight));
    temp_genes = [];
    for i = 1:length(genes_change)
        for j = 1:length(genes_change{i})
            temp_genes = [temp_genes; genes_change{i}(j)];
        end
    end
    genes_change = temp_genes;
    genes_change = unique(genes_change,'stable');

    % fit ellipse to points if = 3 pts available
%     if numSamples == 3
%         binArea = zeros(numBins,1);
%         
%         for iBin = 1:numBins
%             % fit ellipse
%             fprintf('Fitting ellipse to pts: %.2f%%\n',(iBin/numBins*100))
%             
%             % plot large distance pts
%             if ismember(iBin,labelLocs(1:numPts_highlight))
%                 [X, binArea(iBin)] = fitEllipse3d(squeeze(scoreReshape(iBin,1:2,:))');
%                 plot(X(:,1),X(:,2),'r-')
%             end
%         end
%         
%         % fit ellipsoid to points if > 3 pts available
%     elseif numSamples > 3
%         binVol = zeros(numBins,1);
%         
%         for iBin = 1:numBins
%             fprintf('Fitting ellipsoid to pts: %.2f\n',(iBin/numBins*100))
%             
%             % plot large distance pts
%             if ismember(iBin,labelLocs(1:numPts_highlight))
%                 % fit ellipse
%                 [A,C] = MinVolEllipse(squeeze(scoreReshape(iBin,1:3,:)),.01);
%                 
%                 [~,D,~] = svd(A);
%                 a = 1/sqrt(D(1,1));
%                 b = 1/sqrt(D(2,2));
%                 c = 1/sqrt(D(3,3));
%                 binVol(iBin) = (4/3)*pi*a*b*c;
%                 Ellipse_plot(A,C,[.5 0 0]);
%             end
%         end
%     end
    
    %custom legend
    h = zeros(numSamples, 1);
    for iS = 1:numSamples
        h(iS) = plot(NaN,NaN,'.','color',color_opts(iS,:),'markersize',30);
        sampleName{iS} = ['Sample',num2str(iS)];
    end
    legend(h,sampleName);
    


    % format output
    box on
    title(sprintf('Structure-Function Feature Space - %s',dimReducLabel))

    switch dimReduc
        case 'pca'
            xlabel('PC 1'), ylabel('PC 2'), zlabel('PC 3')
        otherwise
            xlabel('Component 1'), ylabel('Component 2'), zlabel('Component 3')
    end
    view(2)
    set(gca,'linewidth',2,'fontsize',15)

elseif strcmp(sfType,'phaseplane')
        %% SF Plot
    clear x y
    numSamples = size(hic,3);
    for iS = 1:numSamples
        [~, pca_score{iS}] = pca((hic(:,:,iS)));
        x(:,iS) = pca_score{iS}(:,1);
        y(:,iS) = log2(feature1D(:,iS)+1);
    end
    idx = sum(log2(feature1D)+1 > 0,2)>0;
    x = x(idx,:);
    y = y(idx,:);

    allDist = zeros(size(x,1),numSamples-1);
    for iS = 1:numSamples-1
        allDist(:,iS) = diag(pdist2([x(:,iS) y(:,iS)], [x(:,iS+1) y(:,iS+1)]));
    end

    numPts_highlight = 10;
    [~,labelLocs] = sort(sum(allDist,2),'descend');


    shape_list = {'o','s','^'};
    colorScale = jet(numSamples);
    figure('position',[50 50 750 550]), hold on
    for iS = 1:numSamples
    %     scatter(x(:,iS),y(:,iS),15,colorScale(iS,:),...
    %         'filled','MarkerFaceAlpha',.25,'MarkerEdgeColor','k','MarkerEdgeAlpha',.1);
        % highlight regions of largest change
    %     scatter(x(labelLocs(1:numPts_highlight),iS), ...
    %              y(labelLocs(1:numPts_highlight),iS), ...
    %              50,colorScale(iS,:),'filled','MarkerEdgeColor','k','Marker',shape_list{iS});

    %     if iS < numSamples
    %     plot([x(labelLocs(1:numPts_highlight),iS),x(labelLocs(1:numPts_highlight),iS+1)]',...
    %         [y(labelLocs(1:numPts_highlight),iS),y(labelLocs(1:numPts_highlight),iS+1)]',...
    %         '-','Color',[.5 .5 .5])
    %     end
    end

    for i = 1:numPts_highlight
            scatter(x(labelLocs(i),1), y(labelLocs(i),1),90,color_opts(i,:),'filled','Marker','o')
            scatter(x(labelLocs(i),2), y(labelLocs(i),2),90,color_opts(i,:),'filled','Marker','s')
            scatter(x(labelLocs(i),3), y(labelLocs(i),3),85,color_opts(i,:),'filled','Marker','^')

        [Xhigh, xArea] =  fitEllipse2d([x(labelLocs(i),:);...
                                        y(labelLocs(i),:)]');
    %      plot(Xhigh(:,1),Xhigh(:,2),'Color',ten_colors(i))

         h = fill(Xhigh(:,1),Xhigh(:,2),color_opts(i,:));
         set(h,'facealpha',.25,'EdgeColor',color_opts(i,:),'LineStyle','-','LineWidth',1.5)
    end

    %custom legend
    % h = zeros(numSamples, 1);
    % for iS = 1:numSamples
    %     h(iS) = plot(NaN,NaN,'.','color',colorScale(iS,:),'markersize',30);
    %     sampleName{iS} = ['Sample',num2str(iS)];
    % end
    % legend(h,sampleName);

    % format output
    box on
    title(sprintf('%s','Phase Portrait'))
    set(gca,'linewidth',2,'fontsize',15,'TickLength',[0 0])
    axis square
    xlabel('Structure (PC1)')
    ylabel('Function')
    set(gca,'XLim',get(gca,'XLim')*1.2,'YLim',get(gca,'YLim')*1.2)
end

end


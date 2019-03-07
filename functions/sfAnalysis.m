function [features,score] = sfAnalysis(hic,rnaSeq,binNames,norm,binInfo,graphWeighted)
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
%
%   Output
%   features:       Concatenated feature array (Nx3 double)
%   score:          Location of bin in low dimensional projection (Nx3 double)
%
%   Scott Ronquist, 1/23/19

%% set default parameters
if nargin<6; graphWeighted=1; end
if nargin<5; binInfo=[]; end
if nargin<4; norm=0; end
if nargin<3; binNames=cellstr(num2str([1:length(rnaSeq)]')); end

%% normalize, remove centromere and extract centrality (depreciated)
if norm
    error(['norm = 1 is currently a depreciated option, normalize Hi-C',...
        'matrices prior to input'])
    hic = norm_hic_bins(hic,binInfo);
end

%% extract centrality
for tp  = 1:size(hic,3)
    G = graph(hic(:,:,tp),'OmitSelfLoops');
    
    % from Sijia, inverse weights for closeness and betweeness
    Ht_root = nthroot(hic(:,:,tp),2);
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

%% add RNA-seq
features = cat(2,features,reshape(rnaSeq,size(rnaSeq,1),1,size(rnaSeq,2)));
Xnorm = FeatureNorm2(features);

% stack time points
XnormStacked = [];
for i = 1:size(Xnorm,3)
    XnormStacked = [XnormStacked;Xnorm(:,:,i)];
end

% PCA
[~,score,~,~,~] = pca(XnormStacked);

%% distance measurement
switch size(Xnorm,3)
    case 1
        
    case 2
        %% CASE 2
        nBins = length(binNames);
        pcaDist = zeros(nBins,1);
        for i = 1:nBins
            fprintf('extracting PCA length: %d\n',i)
            pcaDist(i) = pdist(score([i,nBins+i],1:3));
        end
        
        %plot
        figure, hold on
        for i = 1:nBins
            plot(score(i,1),score(i,2),'bo')
        end
        for i = 1:nBins
            plot(score(i+nBins,1),score(i+nBins,2),'r*')
        end
        
        % select top percentile of genes to label and draw line
        for i = 1:nBins
            if pcaDist(i) > prctile(pcaDist,90)
                plot([score(i,1),score(i+nBins,1)],[score(i,2),score(i+nBins,2)],'k-')
                text(score(i+nBins,1),score(i+nBins,2),binNames{i})
            end
        end
        xlabel('PC1'),ylabel('PC2'),title('S-F features PCA')
        
        %custom legend
        h = zeros(2, 1);
        h(1) = plot(NaN,NaN,'bo');
        h(2) = plot(NaN,NaN,'r*');
        legend(h,'sample1','sample2');
        
    otherwise                   % dimension >= 3
        % calculate ellipse
        nBins = size(binNames,1);
        binVol = zeros(nBins,1);
        temp = [];
        
        % ellipsoid volume
        for i = 1:nBins
            fprintf('extracting bin volume: %d\n',i)
            if size(unique(score(i:nBins:nBins*(size(Xnorm,3)-1)+i,1:3),'rows'),1) > 1
                [A,C] = MinVolEllipse(score(i:nBins:nBins*3+i,1:3)',.03);
                [~,D,~] = svd(A);
                a = 1/sqrt(D(1,1));
                b = 1/sqrt(D(2,2));
                c = 1/sqrt(D(3,3));
                binVol(i) = (4/3)*pi*a*b*c;
                
                % this throwsout ellipsoids that are artificially high (points lie in a plane)
                if D(3,3)/sum(diag(D)) < 1E-10      %D(3,3) < 1e-5
                    binVol(i) = 0;
                end
            end
        end
        
        big_bins = binVol>prctile(binVol,50);
        
        %% plot ellipsoid
        figure, hold on,xlabel('PC1'),ylabel('PC2'),zlabel('PC3'),view(3)
        colormapJet = jet;
        for i = 1:nBins
            fprintf('plotting ellipsoid: %d\n',i)
            if big_bins(i)
                [A ,C] = MinVolEllipse(score(i:nBins:nBins*(size(Xnorm,3)-1)+i,1:3)',.01);
                Ellipse_plot(A,C,colormapJet(ceil((i/nBins)*length(colormapJet)),:));
                text(C(1),C(2),C(3),binNames(i,1),'Color','black','fontsize',15)
            end
        end
        
        %format2
        view(gca,[0 90]);
        set(gca,'fontsize',15,'linewidth',3)
end

end


function [features,score] = sfAnalysis(hic,rnaSeq,binNames,norm,binInfo,graphWeighted)
%sfAnalysis analyzes Hi-C and RNA-seq through centrality and PCA
%   sfAnalysis extracts centrality features from Hi-C, concatonates them
%   with RNA-seq, then determines a low dimensional projection to extract
%   regions that change significantly
%
%   inputs
%   hic: Hi-C matrix, typically normalized (NxNxM double; default: N/A)
%   rnaSeq: rnaseq values, (NxM double; default: N/A)
%   binNames: names associated with Hi-C and RNA-seq bins (string or cell array; default: empty cell)
%   norm: normalize Hi-C flag ([0,1]; default: 0)
%   binInfo: bin chr information necessary for normalization (array; default: [])
%   graphWeighted: use graph weighted centrality analysis ([0,1]; default: 1)
%
%   outputs
%   features: concatonated feature array (Nx3 double)
%   score: location of bin in low dimensional projection (Nx3 double)
%
%   example
%   [B] = function(A)
%
%   Scott Ronquist, 7/28/18

if nargin<6; graphWeighted=1; end
if nargin<5; binInfo=[]; end
if nargin<4; norm=0; end
if nargin<3; binNames=cellstr(num2str([1:length(rnaSeq)]')); end

%% normalize (if necessary), remove centromere and extract centrality
if norm
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
Xnorm_stacked = [];
for i = 1:size(Xnorm,3)
    Xnorm_stacked = [Xnorm_stacked;Xnorm(:,:,i)];
end

% PCA
[~,score,~,~,~] = pca(Xnorm_stacked);

%% distance measurement
switch size(Xnorm,3)
    case 1
        
    case 2
        %% CASE 2
        n_bins = length(binNames);
        pca_dist = zeros(n_bins,1);
        for i = 1:n_bins
            fprintf('extracting PCA length: %d\n',i)
            pca_dist(i) = pdist(score([i,n_bins+i],1:3));
        end
        
        %plot
        figure, hold on
        for i = 1:n_bins
            plot(score(i,1),score(i,2),'bo')
        end
        for i = 1:n_bins
            plot(score(i+n_bins,1),score(i+n_bins,2),'r*')
        end
        
        % select top percentile of genes to label and draw line
        for i = 1:n_bins
            if pca_dist(i) > prctile(pca_dist,90)
                plot([score(i,1),score(i+n_bins,1)],[score(i,2),score(i+n_bins,2)],'k-')
                text(score(i+n_bins,1),score(i+n_bins,2),binNames{i})
            end
        end
        xlabel('PC1'),ylabel('PC2'),title('S-F features PCA')
        
        %custom legend
        h = zeros(2, 1);
        h(1) = plot(NaN,NaN,'bo');
        h(2) = plot(NaN,NaN,'r*');
        legend(h,'sample1','sample2');
        
    case 3
        
    otherwise
        %% OTHERWISE
        % calculate ellipse
        n_bins = size(binNames,1);
        bin_vol = zeros(n_bins,1);
        temp = [];
        
        % ellipsoid volume
        for i = 1:n_bins
            fprintf('extracting bin volume: %d\n',i)
            if size(unique(score(i:n_bins:n_bins*(size(Xnorm,3)-1)+i,1:3),'rows'),1) > 1
                [A,C] = MinVolEllipse(score(i:n_bins:n_bins*3+i,1:3)',.03);
                [~,D,~] = svd(A);
                a = 1/sqrt(D(1,1));
                b = 1/sqrt(D(2,2));
                c = 1/sqrt(D(3,3));
                bin_vol(i) = (4/3)*pi*a*b*c;
                
                % this throwsout ellipsoids that are artificially high (points lie in a plane)
                if D(3,3)/sum(diag(D)) < 1E-10      %D(3,3) < 1e-5
                    bin_vol(i) = 0;
                end
            end
        end
        
        big_bins = bin_vol>prctile(bin_vol,50);
        
        %% plot ellipsoid
        figure, hold on,xlabel('PC1'),ylabel('PC2'),zlabel('PC3'),view(3)
        colormap_jet = jet;
        for i = 1:n_bins
            fprintf('plotting ellipsoid: %d\n',i)
            if big_bins(i)
                [A ,C] = MinVolEllipse(score(i:n_bins:n_bins*(size(Xnorm,3)-1)+i,1:3)',.01);
                Ellipse_plot(A,C,colormap_jet(ceil((i/n_bins)*length(colormap_jet)),:));
                text(C(1),C(2),C(3),binNames(i,1),'Color','black','fontsize',15)
            end
        end
        
        %format2
        view(gca,[0 90]);
        set(gca,'fontsize',15,'linewidth',3)
end

end


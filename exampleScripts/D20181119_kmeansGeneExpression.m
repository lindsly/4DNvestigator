% Goal: create RNA-seq k-means clustering script
%
% Chen, Haiming, et al. "Functional organization of the human 4D Nucleome."
% Proceedings of the National Academy of Sciences 112.26 (2015): 8002-8007.
%
% https://www.mathworks.com/help/bioinfo/examples/gene-expression-profile-analysis.html

%% create sample data
clear
close all

numGenes = 5000;
numTPs = 3;
numK = 4;
TPs = 1:3;
TPLabel = {'G1','S','G2'};
numSubplot = ceil(sqrt(numK));

data = rand(numGenes,numTPs);

%% possible
% probably not for this data, possibilities for filtering are given below
% genevarfilter
% genelowvalfilter
% geneentropyfilter

%% possible normalizations
normType = 'normalize';

switch normType
    case 'log2'
        dataNorm = log2(data+.5)+1;
    case 'normalize'
        dataNorm = normalize(data,2);
    case 'none'
        dataNorm = data;
end

%% k-means
[cidx, ctrs] = kmeans(dataNorm,numK,'dist','corr','rep',5,'disp','final');

%% FIGURE: all lines
figure
for c = 1:numK
    subplot(numSubplot,numSubplot,c);
    plot(TPs,dataNorm((cidx == c),:)');
    xticklabels(TPLabel)
    axis tight
end
suptitle('K-Means Clustering of Profiles');

%% FIGURE: profile trends
figure
for c = 1:numK
    subplot(numSubplot,numSubplot,c);
    plot(TPs,ctrs(c,:)');
    xticklabels(TPLabel)
    axis tight
    axis off
end
suptitle('K-Means Clustering of Profiles');

%% FIGURE: smooth shading and error
figure
for c = 1:numK
    subplot(numSubplot,numSubplot,c);
    stdshade(dataNorm((cidx == c),:),.5,'b',[],.1)
    xticklabels(TPLabel)
    axis tight
end
suptitle('K-Means Smooth trend');

%% FIGURE: clustergram
% get color scale info from rnaSeq
% create a default color map ranging from red to green
cLength = 100;
red = [1 0 0];
green = [0 1 0];
colorScale = [linspace(red(1),green(1),cLength)',...
    linspace(red(2),green(2),cLength)',...
    linspace(red(3),green(3),cLength)'];

figure
for c = 1:numK
    subplot(numSubplot,numSubplot,c);
    imagesc(dataNorm((cidx == c),:))
    colormap(colorScale), colorbar
    xticklabels(TPLabel)
    axis tight
end
suptitle('K-Means Clustergrams');

%% FIGURE: subcluster
for iC = 1:numK
    tempData = dataNorm((cidx == iC),:);
    [cidx2, ~] = kmeans(tempData,numK,'dist','corr','rep',5,'disp','final');
    
    figure
    for c = 1:numK
        subplot(numSubplot,numSubplot,c);
        imagesc(tempData((cidx2 == c),:))
        colormap(colorScale), colorbar
        xticklabels(TPLabel)
        axis tight
    end
    suptitle(sprintf('K-Means sub-Clustergrams cluster:%i',iC));
end


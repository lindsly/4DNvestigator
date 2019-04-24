% Goal: perform Larntz Perlman analysis for 23 x 23 correlation matrix:

%% load data/paths
clear
close all
addpath(genpath('E:\MATLAB\gsfat'))

%% select cell types and Hi-C parameters
%%% PARAMETERS vvv
param.binType = 'BP';
param.binSize = 1E6;
param.norm1d = 'KR';
param.norm3d = 'oe';
param.intraFlag = 1;
param.chr = 14;
param.sampleSelect = {'IMR90','HFFc6','H1-hESC'};
%%% PARAMETERS ^^^

fn = {'http://hicfiles.s3.amazonaws.com/hiseq/rpe1/DarrowHuntley-2015/WT-combined.hic';...
    'https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic';...
    'https://hicfiles.s3.amazonaws.com/hiseq/imr90/in-situ/combined.hic';...
    'https://hicfiles.s3.amazonaws.com/hiseq/hmec/in-situ/combined.hic';...
    'https://hicfiles.s3.amazonaws.com/hiseq/nhek/in-situ/combined.hic';...
    'https://hicfiles.s3.amazonaws.com/hiseq/huvec/in-situ/combined.hic';...
    'https://data.4dnucleome.org/files-processed/4DNFIFLJLIS5/@@download/4DNFIFLJLIS5.hic';...
    'https://data.4dnucleome.org/files-processed/4DNFIOX3BGNE/@@download/4DNFIOX3BGNE.hic'};

cellType = {'RPE1 WT';'GM12878';'IMR90';'HMEC';'NHEK';'HUVEC';'HFFc6';'H1-hESC'};
dataInfo = table(cellType,fn);
dataInfo = sortrows(dataInfo,'cellType','ascend');
dataInfo = dataInfo(ismember(dataInfo.cellType,param.sampleSelect),:);

%% load hic
H = cell(height(dataInfo),1);
hicInfo = {};
for iSample = 1:height(dataInfo)
    disp(iSample)
    [H{iSample},hicInfo{iSample}] = hic2mat('observed','none',dataInfo.fn{iSample},...
        'ALL','ALL',param.binType,1E7,param.intraFlag,1);
end

%% normalize data
% find chr start and end
hChr = zeros(22,22,height(dataInfo));
for iSample = 1:height(dataInfo)
    if iSample == 3
        numChr = 22;
        chrSizes = readtable('hg19.chrom.sizes','filetype','text');
        chrStart = [1;cumsum(ceil(chrSizes{:,2}/param.binSize))+1];
    else
        numChr = 22;
        chrSizes = readtable('hg38.chrom.sizes','filetype','text');
        chrStart = [1;cumsum(ceil(chrSizes{:,2}/param.binSize))+1];
    end
    
    chrStart = chrStart(1:23);
    
    for iChr1 = 1:length(chrStart)-1
        for iChr2 = 1:length(chrStart)-1
            hChr(iChr1,iChr2,iSample) = ...
                nansum(nansum(H{iSample}(chrStart(iChr1):chrStart(iChr1+1)-1,...
                chrStart(iChr2):chrStart(iChr2+1)-1)));
        end
    end
    hChr(:,:,iSample) = hChr(:,:,iSample)/sum(sum(hChr(:,:,iSample)));
end

% correlation
hChrCorr = zeros(size(hChr));
for iSample = 1:height(dataInfo)
    hChrCorr(:,:,iSample) = corr(hChr(:,:,iSample));
end

%% visualize data
% 1Mb
figure
for iSample = 1:height(dataInfo)
    subplot(1,height(dataInfo),iSample)
    imagesc(log(H{iSample})), hold on
    for iChr = 1:22
        plot([chrStart(iChr),chrStart(iChr)],[0 length(H{iSample})],'r-')
        plot([0 length(H{iSample})],[chrStart(iChr),chrStart(iChr)],'r-')
    end
    axis square
end
linkaxes(get(gcf,'children'))

% chr level
figure
for iSample = 1:height(dataInfo)
    subplot(1,height(dataInfo),iSample)
    imagesc(hChr(:,:,iSample))
    title(dataInfo.cellType{iSample})
    axis square
end
linkaxes(get(gcf,'children'))

% chr level - Corr
figure
for iSample = 1:height(dataInfo)
    subplot(1,height(dataInfo),iSample)
    imagesc(hChrCorr(:,:,iSample))
    title(dataInfo.cellType{iSample})
    axis square
end
linkaxes(get(gcf,'children'))

%% LP analyze data
% Larntz-Perlman
alphaParam = .95;
plotFlag = 0;
[H0,M,P,S] = larntzPerlman(hChrCorr,size(hChrCorr,1),alphaParam,plotFlag);

% get L-P regions >= top 5th percentile regions
tempXPrctile = 99;
LPRegions = S > prctile(S(:),tempXPrctile);

% figure
hicCMap = [ones(64,1),[1:-1/63:0]',[1:-1/63:0]'];
climMain = [0 prctile(hChr(:),95)];
figure('units','normalized','position',[.05 .05 .9 .8])
for iSample = 1:height(dataInfo)
    ax = subplot(1,height(dataInfo),iSample);
    imagesc(hChr(:,:,iSample)); axis square, hold on
    title(dataInfo.cellType{iSample})
    
    % colorbar - shift colorbar position to X% of normal height
    colormap(ax,hicCMap); caxis(climMain)
    colorbarScaled(.3)
    
    % add Larntz-Perlman ROIs
    addROICircles(LPRegions)
    
    % add VNGE
    L = diag(sum(hChr(:,:,iSample)))-hChr(:,:,iSample);
    Ln = L*(1/trace(L));
    tempEig = eig(Ln);
    tempEig = eig(hChrCorr(:,:,iSample));tempEig = tempEig/(sum(tempEig));
    vnEntropy = -sum(tempEig.*log(tempEig));
    xlabel(sprintf('VNGE = %.2f',vnEntropy))
end
set(get(gcf,'children'),'linewidth',2,'fontsize',15)
linkaxes(get(gcf,'children'))

%% EXTRA
figure
iSample = 1;
ax = subplot(1,height(dataInfo),iSample);
pt1 = get(ax,{'Position','tightinset','PlotBoxAspectRatio'});
pt1{:}
axis square
pt2 = get(ax,{'Position','tightinset','PlotBoxAspectRatio'});
pt2{:}

imagesc(hChr(:,:,iSample)); axis square, hold on
title(dataInfo.cellType{iSample})

% colorbar - shift colorbar position to X% of normal height
colormap(ax,hicCMap); caxis(climMain)
tempCBar = colorbar;
cBarScale = .3;
tempCBar.Position = [tempCBar.Position(1),...
    tempCBar.Position(2)+tempCBar.Position(4)*(1-cBarScale),...
    tempCBar.Position(3), tempCBar.Position(4)*cBarScale];

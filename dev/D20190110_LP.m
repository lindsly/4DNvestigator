% Goal: test Larntz-Perlman method for comparing covariance matrices and
% create figure for displaying results 
%
%   Reference: Koziol, James A., et al. "A graphical technique for
%   displaying correlation matrices." The American Statistician 51.4
%   (1997): 301-304.

%% Load 
clear
close all

restoredefaultpath
addpath(genpath('.'))

%% Hi-C type example
%%% PARAMETERS vvv
hicParam.binType = 'BP';
hicParam.binSize = 1E6;
hicParam.norm1d = 'KR';
hicParam.norm3d = 'oe';
hicParam.intraFlag = 1;
hicParam.chr = 14;
sampleSelect = {'HFFc6','H1-hESC'};
%%% PARAMETERS ^^^

sampleDataLoc = {'http://hicfiles.s3.amazonaws.com/hiseq/rpe1/DarrowHuntley-2015/WT-combined.hic',...
    'https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic',...
    'https://hicfiles.s3.amazonaws.com/hiseq/imr90/in-situ/combined.hic',...
    'https://hicfiles.s3.amazonaws.com/hiseq/hmec/in-situ/combined.hic',...
    'https://hicfiles.s3.amazonaws.com/hiseq/nhek/in-situ/combined.hic',...
    'https://hicfiles.s3.amazonaws.com/hiseq/huvec/in-situ/combined.hic',...
    'https://data.4dnucleome.org/files-processed/4DNFIFLJLIS5/@@download/4DNFIFLJLIS5.hic',...
    'https://data.4dnucleome.org/files-processed/4DNFIOX3BGNE/@@download/4DNFIOX3BGNE.hic'};

sampleNames = {'RPE1 WT','GM12878','IMR90','HMEC','NHEK','HUVEC','HFFc6','H1-hESC'};

sampleDataLoc = sampleDataLoc(ismember(sampleNames,sampleSelect));
sampleNames = sampleNames(ismember(sampleNames,sampleSelect));

%% load Hi-C
H = [];
for iSample = 1:length(sampleNames)
    tempH = hic2mat(hicParam.norm3d,hicParam.norm1d,sampleDataLoc{iSample},...
        hicParam.chr,hicParam.chr,hicParam.binType,hicParam.binSize,hicParam.intraFlag);
    H = padconcatenation_sr(H,tempH,3);
end

% trim
[HtrimAll,badLocs] = hicTrim(H,1,.5);

% correlation
Hcorr = zeros(size(HtrimAll));
for iSample = 1:length(sampleNames)
    Hcorr(:,:,iSample) = corr(HtrimAll(:,:,iSample));
end

%
figure, 

subplot(1,2,1)
imagesc(Hcorr(:,:,1)), axis square
title('Fibroblast, chr14')
caxis([-1 1])
colorbar

subplot(1,2,2)
imagesc(Hcorr(:,:,2)), axis square
title('ESC, chr14')
caxis([-1 1])
colorbar
set(get(gcf,'children'),'linewidth',2,'fontsize',15)


%% Larntz-Perlman
alphaParam = .95;
plotFlag = 0;
[H0,M,P,S] = larntzPerlman(Hcorr,size(Hcorr,1),alphaParam,plotFlag);

% get L-P regions >= top 5th percentile regions
tempXPrctile = 99;
LPRegions = S > prctile(S(:),tempXPrctile);

% figure
hicCMap = [ones(64,1),[1:-1/63:0]',[1:-1/63:0]'];
climMain = [0 prctile(HtrimAll(:),95)];
figure
for iSample = 1:length(sampleNames)
    ax = subplot(1,length(sampleNames),iSample);
    imagesc(HtrimAll(:,:,iSample)); axis square, hold on
    colormap(ax,hicCMap); colorbar, caxis(climMain)
    title(sampleNames{iSample})
    
    % add Larntz-Perlman ROIs
    addROICircles(LPRegions)
    
    % add VNGE
    L = diag(sum(HtrimAll(:,:,iSample)))-HtrimAll(:,:,iSample);
    Ln = L*(1/trace(L));
    tempEig = eig(Ln);
    tempEig = eig(Hcorr(:,:,iSample));tempEig = tempEig/(sum(tempEig));
    vnEntropy = -sum(tempEig.*log(tempEig));
    xlabel(sprintf('VNGE = %.2f',vnEntropy))
end
set(get(gcf,'children'),'linewidth',2,'fontsize',15)
linkaxes(get(gcf,'children'))

%% simple example
[H0,M,P,S] = larntzPerlman(Hcorr,size(Hcorr,1),alphaParam,plotFlag);
tempXPrctile = 99;
LPRegions = S > prctile(S(:),tempXPrctile);
imagesc(Hcorr(:,:,iSample))
addROICircles(LPRegions,'red')



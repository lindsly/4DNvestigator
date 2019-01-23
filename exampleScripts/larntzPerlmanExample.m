%   This script performs the Larntz-Perlman procedure for 3 distinct Hi-C
%   matrices, as described in "4DNvestigator: a toolbox for the analysis of
%   timeseries Hi-C and RNA-seq data"
%
%   Scott Ronquist, scotronq@umich.edu. 1/21/19

%% Script set-up
clear
close all

%% Hi-C type example
% define Hi-C extraction parameters
%%% PARAMETERS vvv
hicParam.binType = 'BP';
hicParam.binSize = 1E6;
hicParam.norm1d = 'KR';
hicParam.norm3d = 'oe';
hicParam.intraFlag = 1;
hicParam.chr = 14;
sampleSelect = {'IMR90','HFFc6','H1-hESC'};
%%% PARAMETERS ^^^

% Public Hi-C data locations
sampleDataLoc = {'http://hicfiles.s3.amazonaws.com/hiseq/rpe1/DarrowHuntley-2015/WT-combined.hic',...
    'https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic',...
    'https://hicfiles.s3.amazonaws.com/hiseq/imr90/in-situ/combined.hic',...
    'https://hicfiles.s3.amazonaws.com/hiseq/hmec/in-situ/combined.hic',...
    'https://hicfiles.s3.amazonaws.com/hiseq/nhek/in-situ/combined.hic',...
    'https://hicfiles.s3.amazonaws.com/hiseq/huvec/in-situ/combined.hic',...
    'https://data.4dnucleome.org/files-processed/4DNFIFLJLIS5/@@download/4DNFIFLJLIS5.hic',...
    'https://data.4dnucleome.org/files-processed/4DNFIOX3BGNE/@@download/4DNFIOX3BGNE.hic'};
sampleNames = {'RPE1 WT','GM12878','IMR90','HMEC','NHEK','HUVEC','HFFc6','H1-hESC'};

% select samples of interest
sampleDataLoc = sampleDataLoc(ismember(sampleNames,sampleSelect));
sampleNames = sampleNames(ismember(sampleNames,sampleSelect));

%% load Hi-C
H = [];
for iSample = 1:length(sampleNames)
    tempH = hic2mat(hicParam.norm3d,hicParam.norm1d,sampleDataLoc{iSample},...
        hicParam.chr,hicParam.chr,hicParam.binType,hicParam.binSize,hicParam.intraFlag);
    H = padconcatenation_sr(H,tempH,3);
end

%% Process Hi-C
% trim regions with low number of contacts
[HtrimAll,badLocs] = hicTrim(H,1,.5);

% find the correlation matrices
Hcorr = zeros(size(HtrimAll));
for iSample = 1:length(sampleNames)
    Hcorr(:,:,iSample) = corr(HtrimAll(:,:,iSample));
end

%% Larntz-Perlman procedure
% perform the Larntz-Perlman procedure on you correlation matrices
alphaParam = .95;
plotFlag = 0;
[H0,P,S] = larntzPerlman(Hcorr,size(Hcorr,1),alphaParam,plotFlag);

% identify regions in 99th percentile L-P regions
tempXPrctile = 99;
LPRegions = S > prctile(S(:),tempXPrctile);

%% figure
hicCMap = [ones(64,1),[1:-1/63:0]',[1:-1/63:0]'];
climMain = [0 prctile(HtrimAll(:),95)];

figure
for iSample = 1:length(sampleNames)
    
    % plot Hi-C
    ax = subplot(1,length(sampleNames),iSample);
    imagesc(HtrimAll(:,:,iSample)); axis square, hold on
    colormap(ax,hicCMap); colorbar, caxis(climMain)
    title(sampleNames{iSample})
    
    % add circles around Larntz-Perlman ROIs
    addROICircles(LPRegions)
    
    % calculate VNGE
    L = diag(sum(HtrimAll(:,:,iSample)))-HtrimAll(:,:,iSample);
    Ln = L*(1/trace(L));
    tempEig = eig(Ln);
    tempEig = eig(Hcorr(:,:,iSample));tempEig = tempEig/(sum(tempEig));
    vnEntropy = -sum(tempEig.*log(tempEig));
    xlabel(sprintf('VNGE = %.2f',vnEntropy))
end
set(get(gcf,'children'),'linewidth',2,'fontsize',15)
linkaxes(get(gcf,'children'))



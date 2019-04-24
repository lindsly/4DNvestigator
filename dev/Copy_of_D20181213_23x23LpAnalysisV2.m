% Goal: perform Larntz Perlman analysis for 23 x 23 correlation matrix:

%% load data/paths
clear
close all
addpath(genpath('E:\MATLAB\gsfat'))

%% select cell types and Hi-C parameters
%%% PARAMETERS vvv
param.binType = 'BP';
%param.binSize = 1E6;
param.norm1d = 'NONE';%'KR';
param.norm3d = 'observed';%'oe';
param.intraFlag = 1;
param.chr = 14;
param.sampleSelect = {'IMR90','HFFc6','H1-hESC'};
%%% PARAMETERS ^^^

fnPath = {'http://hicfiles.s3.amazonaws.com/hiseq/rpe1/DarrowHuntley-2015/WT-combined.hic';...
    'https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic';...
    'https://hicfiles.s3.amazonaws.com/hiseq/imr90/in-situ/combined.hic';...
    'https://hicfiles.s3.amazonaws.com/hiseq/hmec/in-situ/combined.hic';...
    'https://hicfiles.s3.amazonaws.com/hiseq/nhek/in-situ/combined.hic';...
    'https://hicfiles.s3.amazonaws.com/hiseq/huvec/in-situ/combined.hic';...
    'https://data.4dnucleome.org/files-processed/4DNFIFLJLIS5/@@download/4DNFIFLJLIS5.hic';...
    'https://data.4dnucleome.org/files-processed/4DNFIOX3BGNE/@@download/4DNFIOX3BGNE.hic'};

ct = {'RPE1 WT';'GM12878';'IMR90';'HMEC';'NHEK';'HUVEC';'HFFc6';'H1-hESC'};
dataInfo = table(ct,fnPath);
dataInfo = sortrows(dataInfo,'ct','ascend');
dataInfo = dataInfo(ismember(dataInfo.ct,param.sampleSelect),:);

%% check ref genome - remove non-hg19
fprintf('getting hic header information...\n')
hicHeader = cell(height(dataInfo),1);
remSamp = false(height(dataInfo),1);
for iSample = 1:height(dataInfo)
    [hicHeader{iSample}] = readHicHeader(dataInfo.fnPath{iSample});
    
    % remove chromosome 'All'
    hicHeader{iSample}.Chromosomes(strcmpi(hicHeader{iSample}.Chromosomes.chr,'ALL'),:) = [];
    
    % only keep hg19 hi-c
%     if ismember('refGenome',fieldnames(hicHeader{iSample}))
%         if ~strcmp(hicHeader{iSample}.refGenome,'hg19')
%             remSamp(iSample) = 1;
%         end
%     else
%         remSamp(iSample) = 1;
%     end
end
dataInfo(remSamp,:) = [];
hicHeader(remSamp) = [];

%% Create the waitbar and determine intialization properties
waitBar = waitbar(0,'loading Hi-C...','Name','Loading Data');
set(findall(waitBar),'Units', 'normalized');    % Change the Units Property of the figure and all the children
set(waitBar,'Position', [0.25 0.4 0.5 0.08]);   % Change the size of the figure
totalWait = height(dataInfo);
currentWait = 1;
avTime = 0;

%% load hic
H = cell(height(dataInfo),1);
for iSample = 1:height(dataInfo)
    if iSample==1;tic;end
    
    % update load bar
    waitbar(currentWait/totalWait,waitBar,...
        sprintf('Loading Hi-C. Sample: %s, estimated time: %i secs...',...
        dataInfo.ct{iSample},round(avTime*(totalWait-currentWait))));
    
    % load data
    H{iSample} = hic2mat('observed','none',dataInfo.fnPath{iSample},'ALL','ALL',...
        param.binType,max(hicHeader{iSample}.BasePairdelimitedResolutions),...
        param.intraFlag);
    
    % update load bar variable
    currentWait = currentWait+1;
    if iChr==1 && iSample==1;avTime = toc;end % for estimation of time to load
end

%% normalize data
% find chr start and end
hChr = [];
chrStart = cell(height(dataInfo),1);
for iSample = 1:height(dataInfo)
    % remove M chr
    hicHeader{iSample}.Chromosomes(strncmp(hicHeader{iSample}.Chromosomes.chr,'M',1),:) = [];
    
    % get chr names
    chrName = hicHeader{iSample}.Chromosomes.chr;
    
    % get chr start
    chrStart{iSample} = [1;cumsum(ceil(hicHeader{iSample}.Chromosomes.chrLength/...
        max(hicHeader{iSample}.BasePairdelimitedResolutions)))+1];
    
    % get chr length to normalize for length (in 1Mb)
    chrLengthMat = hicHeader{iSample}.Chromosomes.chrLength/1E6;
    chrLengthMat = chrLengthMat*chrLengthMat';
    
    % construct chr-level hi-c
    hChrTemp = zeros(length(chrName));
    for iChr1 = 1:length(chrStart{iSample})-1
        for iChr2 = 1:length(chrStart{iSample})-1
            hChrTemp(iChr1,iChr2) = ...
                nansum(nansum(H{iSample}(chrStart{iSample}(iChr1):chrStart{iSample}(iChr1+1)-1,...
                chrStart{iSample}(iChr2):chrStart{iSample}(iChr2+1)-1)));
        end
    end
    
    % normalize for chr length
    hChrTemp = hChrTemp./chrLengthMat;
    
    % normalize for total counts
    hChrTemp = (hChrTemp/sum(hChrTemp(:)))*1E3;
    
    hChr = cat(3,hChr,hChrTemp);
end

% correlation
hChrCorr = zeros(size(hChr));
for iSample = 1:height(dataInfo)
    hChrCorr(:,:,iSample) = corr(hChr(:,:,iSample));
end

%% visualize data
% 1Mb resolution
figure
for iSample = 1:height(dataInfo)
    subplot(1,height(dataInfo),iSample)
    imagesc(log(H{iSample})), hold on
    for iChr = 1:22
        plot([chrStart{iSample}(iChr),chrStart{iSample}(iChr)],[0 length(H{iSample})],'r-')
        plot([0 length(H{iSample})],[chrStart{iSample}(iChr),chrStart{iSample}(iChr)],'r-')
    end
    axis square
end
% linkaxes(get(gcf,'children'))

% chr level
figure
for iSample = 1:height(dataInfo)
    subplot(1,height(dataInfo),iSample)
    imagesc(hChr(:,:,iSample))
    title(dataInfo.ct{iSample})
    axis square
end
linkaxes(get(gcf,'children'))

% chr level - Corr
figure
for iSample = 1:height(dataInfo)
    subplot(1,height(dataInfo),iSample)
    imagesc(hChrCorr(:,:,iSample))
    title(dataInfo.ct{iSample})
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
    title(dataInfo.ct{iSample})
    
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
title(dataInfo.ct{iSample})

% colorbar - shift colorbar position to X% of normal height
colormap(ax,hicCMap); caxis(climMain)
tempCBar = colorbar;
cBarScale = .3;
tempCBar.Position = [tempCBar.Position(1),...
    tempCBar.Position(2)+tempCBar.Position(4)*(1-cBarScale),...
    tempCBar.Position(3), tempCBar.Position(4)*cBarScale];

%% Description
% This script provides and overview of the methods described in the
% "4DNvestigator: a toolbox for the analysis of timeseries Hi-C and RNA-seq
% data"
%
% Scott Ronquist, scotronq@umich.edu. 4/29/19

%% Set up
clear
close all

%% Select Data set to Load
% Download publicly available time series Hi-C and RNA-seq datasets
% 
% Time Series Hi-C and RNA-seq data available at the following link:
%   https://drive.google.com/open?id=1lSyU-7I0ME3X70Mt_-HjLMtPc-BMKHxm

if exist('myodData.mat','file')==2 && exist('sampleMyodDataIndexTp-48_8_80.xlsx','file')==2
    error(['Download time series Hi-C and RNA-seq data ',...
        '<a href="https://drive.google.com/open?id=1lSyU-7I0ME3X70Mt_-HjLMtPc-BMKHxm">here</a>'])
end

%% Load data through the 4DNvestigator functions
% 4DNvestigator formatted data is available here:
%   https://drive.google.com/open?id=11XQ0CudxRvM5P6aI57yes2h8XXLhxlRr
%
% Alternatively, if you choose to download the .hic and .genes.results
% files, you can process them through the 4DNvestigator codes

if exist('myodData.mat','file')==2
    load('myodData.mat')
elseif exist('sampleMyodDataIndexTp-48_8_80.xlsx','file')==2
    [dataInfo] = fdnLoadUserInput(indexFile);
    [H] = fdnLoadHic(dataInfo,'single');
    [R] = fdnLoadRnaseq(dataInfo,H);
    
    % save data
    save([dataInfo.path.output,dataInfo.delim,dataInfo.projName,'Data.mat'],...
        'H','R','dataInfo','-v7.3')
else
    error(['Download time series Hi-C and RNA-seq data ',...
        '<a href="https://drive.google.com/open?id=1lSyU-7I0ME3X70Mt_-HjLMtPc-BMKHxm">here</a>'])
end

%% Time series differential expression
% run samples through differential expression analysis, with time series
% gene expression pattern analysis
gseaFlag = 1;
rnaseqPatternFlag = 1;
[R, rnaseqPatternAll] = fdnDiffExpGsaa(dataInfo,R,1,1);

%% 4DN Feature Analyzer Example
% select Regions of Interest
% selecting chromosome 11
chrSelect = 11;
goiH = H.s100kb.oeTrim{chrSelect};
goiR = R.s100kb.tpmMeanTrim{chrSelect};
goi = R.s100kb.geneTrim{chrSelect};

% select dimension reduction method and run through 4DN feature analyzer
dimReduc = 'pca';
[features,score] = sfAnalysis(goiH,goiR,goi,[],[],[],dimReduc);

%% Larntz-Perlman Example
% select Regions of Interest
% selecting region surrounding MYOD1 gene location
goi = 'MYH1';
goiFlank = 20;

% get ROI parameters
goiChr = R.TPM.chr(ismember(R.TPM.geneName,goi),:);
goiLoc = round(mean(find(cell2mat(cellfun(@sum, cellfun(@(x) contains(x,goi),...
    R.s100kb.geneTrim{R.TPM.chr(ismember(R.TPM.geneName,goi),:)},...
    'UniformOutput',false), 'UniformOutput',false)))));
goiLocFlank = max([1, goiLoc-goiFlank]):min([size(H.s100kb.oeTrim{goiChr},1), goiLoc+goiFlank]);

% extract ROI
roiHKr = H.s100kb.krTrim{goiChr}(goiLocFlank,goiLocFlank,:);
roiHOe = H.s100kb.oeTrim{goiChr}(goiLocFlank,goiLocFlank,:);
roiHOeL2 = log2(roiHOe);

% set -Inf values to minimum
tempRoiH = roiHOeL2;
tempRoiH(isinf(tempRoiH)) = NaN;
roiHOeL2(isinf(roiHOeL2)) = nanmin(tempRoiH(:));

% Calculate the correlation matrices
Hcorr = zeros(size(roiHOeL2));
for iSample = 1:size(roiHOeL2,3)
    Hcorr(:,:,iSample) = corr(roiHOeL2(:,:,iSample));
end

% Perform the Larntz-Perlman procedure on these correlation matrices
alphaParam = .95;
plotFlag = 0;
[H0,P,S] = larntzPerlman(Hcorr,size(Hcorr,1),alphaParam,plotFlag);

% Identify regions in 99th percentile L-P regions
tempXPrctile = 99;
LPRegions = S > prctile(S(:),tempXPrctile);

% create blue-white-red colormap based on min/max
roiHLimAbs = min(abs([nanmin(tempRoiH(:)) nanmax(tempRoiH(:))]));
roiHLim = [-roiHLimAbs roiHLimAbs];
hicCMap = 1-((1-redblue(100))*.7);

% Figure
figure('position',[100 500 1310 340])
for iSample = 1:size(roiHOeL2,3)
    
    % Plot Hi-C
    ax = subplot(1,size(roiHOeL2,3),iSample);
    imagesc(roiHOeL2(:,:,iSample)), axis square
    colormap(ax, hicCMap)
    colorbar
    caxis(roiHLim)
    
    if iSample==1; ylabel('log_2(O/E Hi-C)'); end
    title(sprintf('%ih', dataInfo.sampleInfo.timePoint(iSample)))
    
    % Add circles around Larntz-Perlman ROIs
    addROICircles(LPRegions,'green',1.5,4)
end
set(get(gcf,'children'),'linewidth',2,'fontsize',15)
linkaxes(get(gcf,'children'))




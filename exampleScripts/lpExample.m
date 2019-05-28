%% Larntz-Perlman (LP) method example
% This example shows how the Larntz-Perlman (LP) method can be used to
% determine where Hi-C matrices are different between samples
%
%   link to paper: (In preparation)
%
%   Version 1.0 (5/24/19)
%   Written by: Scott Ronquist
%   Contact:    scotronq@umich.edu
%   Created:    5/24/19
%
%   Revision History:
%   v1.0 (5/24/19)
%   * lpExample.m created

%% Load Data
% 4DNvestigator formatted Hi-C and RNA-seq data is available here:
%   https://drive.google.com/open?id=11XQ0CudxRvM5P6aI57yes2h8XXLhxlRr
clear
close all

load('myodData.mat')

%% Pre-processing Data
% select Regions of Interest: region surrounding "roi"
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

%% Larntz-Perlman procedure
alphaParam = .05;
plotFlag = 0;
[H0,P,S,pval] = larntzPerlman(Hcorr,size(Hcorr,1),alphaParam,plotFlag);

% Identify regions in 99th percentile L-P regions
tempXPrctile = 99;
LPRegions = S > prctile(S(:),tempXPrctile);

% create blue-white-red colormap based on min/max
roiHLimAbs = min(abs([nanmin(tempRoiH(:)) nanmax(tempRoiH(:))]));
roiHLim = [-roiHLimAbs roiHLimAbs];
hicCMap = 1-((1-redblue(100))*.7);

%% Figure
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

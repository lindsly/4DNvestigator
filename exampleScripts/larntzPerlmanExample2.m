%   This script performs the Larntz-Perlman procedure for 3 distinct Hi-C
%   matrices, in specific regions of interest, as described in
%   "4DNvestigator: a toolbox for the analysis of timeseries Hi-C and
%   RNA-seq data"
%
%   Version 1.0 (4/21/19)
%   Written by: Scott Ronquist
%   Contact: scotronq@umich.edu
%   Contributors:
%   Created: 4/21/19
%   Revision History:
%   v1.0 (4/21/19)
%   * Created script

%% set up
clear
close all

%% load data
if ~isfile('./data/myodTsData.mat')
    [dataInfo] = fdnLoadUserInput('myodDataIndex.xlsx','myod','.');
    [H] = fdnLoadHic(dataInfo);
    [R] = fdnLoadRnaseq(dataInfo,H);
    
    save('./data/myodTsData','H','R','dataInfo','-v7.3')
else
    load('./data/myodTsData.mat')
end

%% set-up for L-P analysis
% select ROIs
goi = 'MYOD1';
goiFlank = 50;
goiChr = R.TPM.chr(ismember(R.TPM.geneName,goi),:);
goiLoc = find(cell2mat(cellfun(@sum, cellfun(@(x) contains(x,goi),...
    R.s100kb.geneTrim{R.TPM.chr(ismember(R.TPM.geneName,goi),:)},...
    'UniformOutput',false), 'UniformOutput',false)));
goiLocFlank = max([1, goiLoc-goiFlank]):min([size(H.s100kb.oeTrim{goiChr},1), goiLoc+goiFlank]);
roiH = H.s100kb.oeTrim{goiChr}(goiLocFlank,goiLocFlank,:);

% Calculate the correlation matrices
Hcorr = zeros(size(roiH));
for iSample = 1:size(roiH,3)
    Hcorr(:,:,iSample) = corr(roiH(:,:,iSample));
end

% figure check
figure
roiHLim = [min(log(roiH(:))) max(log(roiH(:)))];
hicCMap = [ones(64,1),[1:-1/63:0]',[1:-1/63:0]'];
for iSample = 1:size(roiH,3)
    ax = subplot(2,size(roiH,3),iSample);
    imagesc(log(roiH(:,:,iSample))), axis square
    colormap(ax, hicCMap), caxis(roiHLim)
    title(sprintf('log_2(Hi-C), sample %i', iSample))
    
    subplot(2,size(roiH,3),iSample+size(roiH,3))
    imagesc(Hcorr(:,:,iSample)), axis square
    caxis([-1 1])
    title(sprintf('corr(Hi-C), sample %i', iSample))
end
set(get(gcf,'children'),'linewidth',2,'fontsize',15)

%% Larntz-Perlman procedure
% Perform the Larntz-Perlman procedure on your correlation matrices
alphaParam = .95;
plotFlag = 0;
[H0,P,S] = larntzPerlman(Hcorr,size(Hcorr,1),alphaParam,plotFlag);

% Identify regions in 99th percentile L-P regions
tempXPrctile = 99;
LPRegions = S > prctile(S(:),tempXPrctile);

%% Figure
figure('position',[100 500 1310 340])
for iSample = 1:size(roiH,3)
    
    % Plot Hi-C
    ax = subplot(1,size(roiH,3),iSample);
    imagesc(log(roiH(:,:,iSample))), axis square
    colormap(ax, hicCMap), caxis(roiHLim)
    title(sprintf('log_2(Hi-C), sample %i', iSample))
    
    % Add circles around Larntz-Perlman ROIs
    addROICircles(LPRegions)
    
    % Calculate VNGE
    vneOut = vne(Hcorr(:,:,iSample));
    %xlabel(sprintf('VNGE = %.2f',vneOut))
end
set(get(gcf,'children'),'linewidth',2,'fontsize',15)
linkaxes(get(gcf,'children'))



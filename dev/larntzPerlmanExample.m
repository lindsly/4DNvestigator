%   This script performs the Larntz-Perlman procedure for 3 distinct Hi-C
%   matrices, as described in "4DNvestigator: a toolbox for the analysis of
%   timeseries Hi-C and RNA-seq data"
%
%   Version 1.1 (3/7/19)
%   Written by: Scott Ronquist
%   Contact: scotronq@umich.edu
%   Contributors:
%   Created: 1/21/19
%   Revision History:
%   v1.1 (3/7/19)
%   * Updated script with dialog box for sample selection

%% Script set-up
clear
close all

%% load data
% Public Hi-C data locations
sampleDataLoc = {'http://hicfiles.s3.amazonaws.com/hiseq/rpe1/DarrowHuntley-2015/WT-combined.hic',...
    'https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic',...
    'https://hicfiles.s3.amazonaws.com/hiseq/imr90/in-situ/combined.hic',...
    'https://hicfiles.s3.amazonaws.com/hiseq/hmec/in-situ/combined.hic',...
    'https://hicfiles.s3.amazonaws.com/hiseq/nhek/in-situ/combined.hic',...
    'https://hicfiles.s3.amazonaws.com/hiseq/huvec/in-situ/combined.hic',...
    'https://data.4dnucleome.org/files-processed/4DNFIFLJLIS5/@@download/4DNFIFLJLIS5.hic',...
    'https://data.4dnucleome.org/files-processed/4DNFIOX3BGNE/@@download/4DNFIOX3BGNE.hic',...
    'https://data.4dnucleome.org/files-processed/4DNFIZ4F74QR/@@download/4DNFIZ4F74QR.hic'};
sampleNames = {'RPE1 WT','GM12878','IMR90','HMEC','NHEK','HUVEC','HFFc6','H1-hESC','HFF-hTERT'};

%% User input for cell type select
%%%% SETTING DIALOG OPTIONS
title = 'Input Larntz-Perlman Analysis parameters';

% Options.WindowStyle = 'modal';
options.Resize = 'on';
options.Interpreter = 'tex';
options.CancelButton = 'on';
options.ApplyButton = 'on';
options.ButtonNames = {'Continue','Cancel'}; %<- default names, included here just for illustration
options.Dim = 4; % Horizontal dimension in fields

prompt = {};
formats = {};
defAns = struct([]);

% Select cell type sample
prompt(1,:) = {'Select samples for LP analysis:','sample',[]};
formats(1,1).type = 'list';
formats(1,1).style = 'listbox';
formats(1,1).format = 'text'; % Answer will give value shown in items, disable to get integer
formats(1,1).items = sampleNames;
formats(1,1).limits = [0 4]; % multi-select
formats(1,1).size = [140 80];
defAns(1).sample = {'HFFc6','H1-hESC','HFF-hTERT'};

% Select chromosome
prompt(2,:) = {'Select chromosome:','chr',[]};
formats(2,1).type   = 'list';
formats(2,1).style  = 'popupmenu';
formats(2,1).items  = 1:22;

% get user inputs parameters
[answer,Cancelled] = inputsdlg(prompt,title,formats,defAns,options);
hicParam.chr = answer.chr;
sampleSelect = answer.sample;

%% Hi-C type example
% Define Hi-C extraction parameters
%%% PARAMETERS vvv
hicParam.binType = 'BP';
hicParam.binSize = 1E6;
hicParam.norm1d = 'KR';
hicParam.norm3d = 'oe';
hicParam.intraFlag = 1;
%%% PARAMETERS ^^^

% Select samples of interest
sampleDataLoc = sampleDataLoc(ismember(sampleNames,sampleSelect));
sampleNames = sampleNames(ismember(sampleNames,sampleSelect));

%% Load Hi-C
H = [];
for iSample = 1:length(sampleNames)
    tempH = hic2mat(hicParam.norm3d,hicParam.norm1d,sampleDataLoc{iSample},...
        hicParam.chr,hicParam.chr,hicParam.binType,hicParam.binSize,hicParam.intraFlag);
    H = padconcatenation_sr(H,tempH,3);
end

%% Process Hi-C
% Trim regions with a low number of contacts
[HtrimAll,badLocs] = hicTrim(H,1,.5);

% Calculate the correlation matrices
Hcorr = zeros(size(HtrimAll));
for iSample = 1:length(sampleNames)
    Hcorr(:,:,iSample) = corr(HtrimAll(:,:,iSample));
end

%% Larntz-Perlman procedure
% Perform the Larntz-Perlman procedure on your correlation matrices
alphaParam = .95;
plotFlag = 0;
[H0,P,S] = larntzPerlman(Hcorr,size(Hcorr,1),alphaParam,plotFlag);

% Identify regions in 99th percentile L-P regions
tempXPrctile = 99;
LPRegions = S > prctile(S(:),tempXPrctile);

%% Figure
hicCMap = [ones(64,1),[1:-1/63:0]',[1:-1/63:0]'];
climMain = [0 prctile(HtrimAll(:),95)];

figure('position',[100 500 1310 340])
for iSample = 1:length(sampleNames)
    
    % Plot Hi-C
    ax = subplot(1,length(sampleNames),iSample);
    imagesc(HtrimAll(:,:,iSample)); axis square, hold on
    colormap(ax,hicCMap); colorbar, caxis(climMain)
    title(sampleNames{iSample})
    
    % Add circles around Larntz-Perlman ROIs
    addROICircles(LPRegions)
    
    % Calculate VNGE
    L = diag(sum(HtrimAll(:,:,iSample)))-HtrimAll(:,:,iSample);
    Ln = L*(1/trace(L));
    tempEig = eig(Ln);
    tempEig = eig(Hcorr(:,:,iSample));tempEig = tempEig/(sum(tempEig));
    vnEntropy = -sum(tempEig.*log(tempEig));
    xlabel(sprintf('VNGE = %.2f',vnEntropy))
end
set(get(gcf,'children'),'linewidth',2,'fontsize',15)
linkaxes(get(gcf,'children'))



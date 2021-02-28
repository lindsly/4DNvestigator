%% LP method benchmarking
% This example shows how the LP method compares with alternate Hi-C
% comparison methods
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
%   * lpBenchmarkExample.m created

%% Set up
clear
close all

% Add 4DNvestigator tools to path
filepath = mfilename('fullpath');
fdnPath = filepath(1:strfind(filepath,'4DNvestigator')+12);
addpath(genpath(fdnPath))

%% Load data samples
% samples and file locations
[dataInfo] = fdnLoadUserInput('benchmarkHicDataIndex.xlsx','lpBenchmark',...
    fullfile('.',filesep,'lpBenchmarkOutput'));

% Hi-C header check
hicHeader = cell(height(dataInfo.sampleInfo),1);
for iSample = 1:height(dataInfo.sampleInfo)
    hicHeader{iSample} = readHicHeader(dataInfo.sampleInfo.path{iSample},dataInfo.sampleInfo.refGenome{iSample});
end

% genome information - keep HGNC synmbols
dataInfo.biomart = readtable('mart_export_ensembl_hg37_info.txt');
[C,IA,IC] = unique(dataInfo.biomart.HGNCSymbol);
dataInfo.biomart = dataInfo.biomart(IA,:);

%% Select ROI parameters
%%%%%%%%%%%%%
roiGene = 'MYOD1';          % GOI to simulate from
roiSamples = 1;             % select sample from "dataInfo.sampleInfo" to construct simulated data from
roiBinsize = 10E3;          % 10E3,50E3
roiBinFlank = 20;           % # of bins flanking GOI
simulationType = 'compartment';    % compartment, loop
%%%%%%%%%%%%%

% get roiGene TSS bin
roiChrTss = [dataInfo.biomart.Chromosome_scaffoldName(ismember(dataInfo.biomart.HGNCSymbol,roiGene)),...
    dataInfo.biomart.GeneStart_bp_(ismember(dataInfo.biomart.HGNCSymbol,roiGene))];

%% Extract Hi-C
roiHKr = [];
roiHOe = [];
for iSample = 1:length(roiSamples)
    fprintf('loading Sample: %i...\n',iSample)
    
    % get ROI locs
    roiLoc = sprintf('%s:%i:%i',roiChrTss{1},...
        roiChrTss{2}-(roiBinsize*roiBinFlank),...
        roiChrTss{2}+(roiBinsize*roiBinFlank));
    
    % extract from Hi-C - KR
    temp = juicerToolsDump('observed','kr',...
        dataInfo.sampleInfo.path{roiSamples(iSample)},...
        roiLoc,roiLoc,'BP',roiBinsize);
    temp(:,1:2) = temp(:,1:2)-min(temp(:,1))+1;
    temp2 = zeros(max(temp(:,2)));
    temp2(sub2ind([max(temp(:,2)) max(temp(:,2))], temp(:,1), temp(:,2))) = temp(:,3);
    temp2 = temp2 + triu(temp2,1)';
    roiHKr = padconcatenation_sr(roiHKr,temp2,3);
    
    % extract from Hi-C - OE
    temp = juicerToolsDump('oe','kr',...
        dataInfo.sampleInfo.path{roiSamples(iSample)},...
        roiLoc,roiLoc,'BP',roiBinsize);
    temp(:,1:2) = temp(:,1:2)-min(temp(:,1))+1;
    temp2 = zeros(max(temp(:,2)));
    temp2(sub2ind([max(temp(:,2)) max(temp(:,2))], temp(:,1), temp(:,2))) = temp(:,3);
    temp2 = temp2 + triu(temp2,1)';
    roiHOe = padconcatenation_sr(roiHOe,temp2,3);
end
roiHOeL2 = log2(roiHOe);

%% Create Simulation data matrices
sampSelect = 1;
switch simulationType
    case 'loop'
        % loop data simulation
        % add "counts" to matrix in specific location
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        countLoc = sort([10 30],'ascend');%sort([5 10],'ascend');           % location of additional counts
        numMats = 10;           % number of simulated matrices, incremental up to max
        totPcrt = 6;%4;         % times above mean that is added
        rng(1)                  % random number generator
        gaussAddFlag = 1;       % make additiona gaussian
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % get Expected matrix, percentages to add, and find average count for diag
        roiHAddE = roiHKr(:,:,sampSelect)./roiHOe(:,:,sampSelect);
        prct2add = linspace(0, totPcrt, numMats);
        meanCount = mean(diag(roiHKr(:,:,sampSelect),diff(countLoc)));
        N = size(roiHAddE,1);
        
        % create synthetic matrices
        roiHAddKr = zeros([size(roiHAddE), numMats]);
        roiHAddOe = zeros([size(roiHAddE), numMats]);
        roiHAddOeL2 = zeros([size(roiHAddE), numMats]);
        
        gaussAddMat = zeros(size(roiHAddE)); % figure, imagesc(gaussAddMat)
        gaussAddMat(countLoc(1), countLoc(2)) = 1; gaussAddMat(countLoc(2), countLoc(1)) = 1;
        gaussAddMat = imgaussfilt(gaussAddMat,1);
        
        for iMat = 1:numMats
            % add counts to location in KR
            roiHAddKr(:,:,iMat) = roiHKr(:,:,sampSelect);
            if gaussAddFlag
                roiHAddKr(:,:,iMat) = roiHAddKr(:,:,iMat) + gaussAddMat.*(meanCount*prct2add(iMat)*2);
            else
                roiHAddKr(countLoc(1),countLoc(2),iMat) = ...
                    roiHAddKr(countLoc(1),countLoc(2),iMat)+(meanCount*prct2add(iMat));
                roiHAddKr(countLoc(2),countLoc(1),iMat) = ...
                    roiHAddKr(countLoc(2),countLoc(1),iMat)+(meanCount*prct2add(iMat));
            end
            
            % get new O/E
            roiHAddOe(:,:,iMat) = roiHAddKr(:,:,iMat)./roiHAddE;
            roiHAddOeL2(:,:,iMat) = log2(roiHAddOe(:,:,iMat));
        end
        
    case 'compartment'
        % Create Simulated Compartment Change Data
        % change A/B structure "counts" to matrix in specific location
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        changeLoc = 16;
        numMats = 10;
        rng(1)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % get Expected matrix, percentages to add, and find average count for diag
        roiHAddE = roiHKr(:,:,sampSelect)./roiHOe(:,:,sampSelect);
        corrVals = linspace(1, -.5, numMats+1);corrVals(1) = [];
        N = size(roiHAddE,1);
        
        % create synthetic matrices
        roiHAddKr = zeros([size(roiHAddE), numMats]);
        roiHAddOe = zeros([size(roiHAddE), numMats]);
        roiHAddOeL2 = zeros([size(roiHAddE), numMats]);
        tempRho = zeros(numMats,1);
        for iMat = 1:numMats
            
            % change log2 O/E correlation at changeLoc
            roiHAddOeL2(:,:,iMat) = roiHOeL2(:,:,sampSelect);
            
            % create random correlated data
            % https://stackoverflow.com/questions/16713469/generating-two-correlated-random-vectors
            M = [roiHAddOeL2(:,changeLoc,iMat) std(roiHAddOeL2(:,changeLoc,iMat))*randn(N,1)];
            L = chol([1 corrVals(iMat); corrVals(iMat) 1]);
            M = M*L;
            
            % smooth signal (creates data more similar to real Hi-C)
            M(:,2) = movmean(M(:,2),3);
            
            % add counts to synthetic data
            roiHAddOeL2(:,changeLoc,iMat) = M(:,2);
            roiHAddOeL2(changeLoc,:,iMat) = M(:,2)';
            
            % get correlation after smoothing and noise
            tempRho(iMat) = corr(roiHAddOeL2(:,changeLoc,iMat), roiHOeL2(:,changeLoc,sampSelect));
            
            % recreate KR and O/E
            roiHAddOe(:,:,iMat) = 2.^roiHAddOeL2(:,:,iMat);
            roiHAddKr(:,:,iMat) = roiHAddOe(:,:,iMat).*roiHAddE;
        end
        
        % add a noise matrix
        tempN = triu(randn(N,N))*.05; tempN = tempN + triu(tempN,1)';
        roiHAddOeL2 = cat(3,roiHAddOeL2,roiHOeL2(:,:,sampSelect) + tempN);
        roiHAddOe = cat(3,roiHAddOe,2.^roiHAddOeL2(:,:,end));
        roiHAddKr = cat(3,roiHAddKr,roiHAddOe(:,:,end).*roiHAddE);
        numMats = numMats+1;
        tempRho = [tempRho;corr(roiHAddOeL2(:,changeLoc,end),roiHOeL2(:,changeLoc,sampSelect))];
    otherwise
        error('please select valid "simulationType" variable: "loop" or "compartment"')
end

%% Figure to show simulated data matrices
hicCMap = 1-((1-redblue(100))*.7);
figure,
subplot(4,3,1)
imagesc(roiHOeL2(:,:,sampSelect)), axis square
colormap(hicCMap), caxis([-2 2]), colorbar
title('input Hi-C matrix')
ylabel('log_2(Hi-C O/E)')

for iMat = 1:numMats
    subplot(4,3,iMat+1)
    imagesc(roiHAddOeL2(:,:,iMat)), axis square
    colormap(hicCMap), caxis([-2 2])
    switch simulationType
        case 'loop'
            title(sprintf('Mat:%i, counts added:%i',iMat,...
                round(meanCount*prct2add(iMat))))
        case 'compartment'
            title(sprintf('Mat:%i, corr:%.2f',iMat,tempRho(iMat)))
    end
end

linkaxes(get(gcf,'children'))
set(get(gcf,'children'),'linewidth',2)

%% Test Multiple Hi-C comparison methods
% test diff Hi-C measures
PvalLP = zeros(numMats,1);
PvalHicSpector = zeros(numMats,1);
PvalHicRep = zeros(numMats,1);
PvalSelfish = zeros(numMats,1);
Ps = cell(numMats,1);
for iMat = 1:numMats
    % get test matrices - original concatenated with test
    tempRoiHOe = cat(3,roiHOe(:,:,sampSelect),roiHAddOe(:,:,iMat));
    tempRoiHKr = cat(3,roiHKr(:,:,sampSelect),roiHAddKr(:,:,iMat));
    roiHOeL2 = log2(tempRoiHOe);
    
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
    alphaParam = .05;
    plotFlag = 0;
    [H0,P,S,PvalLP(iMat)] = larntzPerlman(Hcorr,size(Hcorr,1),alphaParam,plotFlag);
    
    % HicSpector
    [~,PvalHicSpector(iMat)] = hicSpector(tempRoiHKr(:,:,1),tempRoiHKr(:,:,2),20);
    
    % HiCRep
    [PvalHicRep(iMat)] = hicRep(tempRoiHKr(:,:,1:2),roiBinsize);
    
    % SELFISH
    THRESHOLD = .1; RESOLUTION = 1; INTERVAL = [1 40];
    temp1 = zeros(nchoosek(N,2)+N,3); temp2 = zeros(nchoosek(N,2)+N,3);count = 1;
    for i=1:N
        for j=i:N
            temp1(count,:) = [i,j,tempRoiHKr(i,j,1)];
            temp2(count,:) = [i,j,tempRoiHKr(i,j,2)];
            count = count+1;
        end
    end
    [X,Y,Ps{iMat}] = SELFISH_DCI(temp1,temp2,[],[],THRESHOLD,RESOLUTION,INTERVAL);
    X = X+1; Y = Y+1;
    if isempty(Ps{iMat})
        PvalSelfish(iMat) = 1;
    else
        PvalSelfish(iMat) = min(Ps{iMat});
    end
    
    % figure to show location of SELFISH DCI detection
    if 1==0
        figure
        subplot(2,2,1), imagesc(log(tempRoiHKr(:,:,1))), axis square, hold on, plot(X,Y,'r*')
        subplot(2,2,2), imagesc(log(tempRoiHKr(:,:,2))), axis square, hold on, plot(X,Y,'r*')
        tempP = ones(N); tempP(sub2ind([N N], X, Y)) = Ps{iMat};
        subplot(2,2,3), imagesc(tempP), axis square, caxis([0 .1])
        subplot(2,2,4), imagesc(S), axis square
        linkaxes(get(gcf, 'children'))
    end
end

% plot performance as matrices diverge
figure, plot(1:numMats, [PvalLP PvalHicSpector PvalHicRep PvalSelfish],'*-','linewidth',2)
xlim([1 numMats]), ylim([0 1])
xticks(1:numMats)
xlabel('Matrix Comparison'), ylabel('Difference Measure')
legend({'L-P','HiC-spector','HiCRep','SELFISH'})
set(gca,'linewidth',2,'fontsize',15)

switch simulationType
    case 'loop'
        xticklabels(round(meanCount*prct2add))
        xlabel('Counts added')
    case 'compartment'
        xticklabels(round(tempRho, 2)), xtickangle(45)
        xlabel('Correlation change')
end
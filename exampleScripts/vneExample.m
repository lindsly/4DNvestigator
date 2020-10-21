function [] = vneExample(Folder_Data, Folder_Result, chrSelect, bpFrag, binSize)
    %% Quantifying Structural Order through Network Entropy (VNE)
    % This example shows how von Neumann Entropy (VNE) can be used to quantify
    % chromatin structural order in Hi-C data.
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
    %   * vneExample.m created

    %% Load Data
    % clear
    % close all
    % 
    % % Add 4DNvestigator tools to path
    % filepath = mfilename('fullpath');
    % fdnPath = filepath(1:strfind(filepath,'4DNvestigator')+12);
    % addpath(genpath(fdnPath))

    % paths to processed HFFc6 and H1-hESC data
    % fn = {'https://data.4dnucleome.org/files-processed/4DNFIFLJLIS5/@@download/4DNFIFLJLIS5.hic',...
    %     'https://data.4dnucleome.org/files-processed/4DNFIOX3BGNE/@@download/4DNFIOX3BGNE.hic'};

    % fn = {'https://drive.google.com/file/d/1ezQtUFWQMKAxry4web5Zm7mMilEUaw_6/view?usp=sharing',...
    %     'https://drive.google.com/file/d/1YsEgqRv8deXRxAujj_i8R4jT9TqMbedq/view?usp=sharing'};

    % TEMPORARY
    fn = {'\\172.17.109.24\internal_4dn\projects\4DNvestigator_data\4DNFIFLJLIS5.hic',...
          '\\172.17.109.24\internal_4dn\projects\4DNvestigator_data\4DNFIOX3BGNE.hic'};
    % Hi-C matrix parameters

    %% Set default 4DN Feature Analyzer parameters
    if ~exist('chrSelect','var')||isempty(chrSelect);chrSelect = 14;end
    if ~exist('bpFrag','var')||isempty(bpFrag);bpFrag='BP';end
    if ~exist('binSize','var')||isempty(binSize);binSize=1E5;end

    % load Hi-C data
    hHFFc6 = hic2mat('oe','kr',fn{1},chrSelect,chrSelect,bpFrag,binSize,1,0);
    hESC = hic2mat('oe','kr',fn{2},chrSelect,chrSelect,bpFrag,binSize,1,0);

    %% Process data
    % remove disconnected nodes (genomic regions with a low number of aligned
    % reads)
    [HTrim,badLocs] = hicTrim(cat(3,hHFFc6,hESC),2,.1);

    % compute the log2, change -inf values to minimum
    HTrim = log2(HTrim);
    for iA = 1:length(fn)
        tempH = HTrim(:,:,iA);
        tempH(tempH==-inf) = min(tempH(isfinite(tempH)));
        HTrim(:,:,iA) = tempH;
    end

    %% VNE Computation
    % set VNE parameters
    preProcess = 'corr';

    % calculate VNE
    vnEntropy = hicVnEntropy(HTrim,[],[],preProcess);

    %% Visualize Matrices with VNE
    % figure
    figure('position',[50 50 1300 500])
    hicCMap = 1-((1-redblue(100))*.7);

    subplot(1,2,1)
    imagesc(HTrim(:,:,1)), axis square
    title(sprintf('HFFc6, VNE: %.2f',vnEntropy(1)))
    colormap(hicCMap), caxis([-2 2]), colorbar
    ylabel('log_2(O/E)')

    subplot(1,2,2)
    imagesc(HTrim(:,:,2)), axis square
    title(sprintf('H1-hESC, VNE: %.2f',vnEntropy(2)))
    colormap(hicCMap), caxis([-2 2])

    set(get(gcf,'children'),'linewidth',2,'fontsize',20)
    linkaxes(get(gcf,'children'))

    %% Visualize Matrices with VNE - corr
    % figure
    figure('position',[50 50 1300 500])
    hicCMap = 1-((1-redblue(100))*.7);

    subplot(1,2,2)
    imagesc(corr(HTrim(:,:,1))), axis square
    caxis([-1 1])

    subplot(1,2,1)
    imagesc(corr(HTrim(:,:,2))), axis square
    caxis([-1 1]), colorbar
    ylabel('corr(log_2(O/E))')

    set(get(gcf,'children'),'linewidth',2,'fontsize',20)
    linkaxes(get(gcf,'children'))

    %% Visualize Matrices with VNE
    % figure
    figure('position',[50 50 1300 500])
    hicCMap = 1-((1-redblue(100))*.7);

    subplot(1,2,2)
    imagesc(HTrim(:,:,1)), axis square
    colormap(hicCMap), caxis([-2 2])

    subplot(1,2,1)
    imagesc(HTrim(:,:,2)), axis square
    colormap(hicCMap), caxis([-2 2]), colorbar
    ylabel('log_2(O/E)')

    set(get(gcf,'children'),'linewidth',2,'fontsize',20)
    linkaxes(get(gcf,'children'))

end
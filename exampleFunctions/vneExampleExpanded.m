function [] = vneExampleExpanded(Data_Loc, Folder_Result, chrSelect, bpFrag, binSize)
    % Quantifying Structural Order through Network Entropy (VNE)
    % This example shows how von Neumann Entropy (VNE) can be used to quantify
    % chromatin structural order in Hi-C data. This example includes addtional
    % samples and cell types to demonstrate the robustness of this measure.
    %
    %   link to paper: (In preparation)
    %
    %   Version 1.0 (8/13/19)
    %   Written by: Scott Ronquist
    %   Contact:    scotronq@umich.edu
    %   Created:    8/13/19
    %
    %   Revision History:
    %   v1.0 (8/13/19)
    %   * vneExampleExpanded.m created

    %% Load Data
    % clear
    % close all

    % Add 4DNvestigator tools to path
    % filepath = mfilename('fullpath');
    % fdnPath = filepath(1:strfind(filepath,'4DNvestigator')+12);
    % addpath(genpath(fdnPath))

    % Paths to processed cell type data
    % Data_Loc = { 'https://data.4dnucleome.org/files-processed/4DNFIFLJLIS5/@@download/4DNFIFLJLIS5.hic';...
    %             'https://data.4dnucleome.org/files-processed/4DNFIOX3BGNE/@@download/4DNFIOX3BGNE.hic';...
    %             'https://hicfiles.s3.amazonaws.com/hiseq/imr90/in-situ/combined.hic';...
    %             'https://hicfiles.s3.amazonaws.com/hiseq/huvec/in-situ/combined.hic';...
    %             'https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic';...
    %             'http://hicfiles.s3.amazonaws.com/hiseq/rpe1/DarrowHuntley-2015/WT-combined.hic'};
    % TEMPORARY
%     Data_Loc = { '\\172.17.109.24\internal_4dn\projects\4DNvestigator_data\4DNFIFLJLIS5.hic';...
%                 '\\172.17.109.24\internal_4dn\projects\4DNvestigator_data\4DNFIOX3BGNE.hic';...
%                 'https://hicfiles.s3.amazonaws.com/hiseq/imr90/in-situ/combined.hic';...
%                 'https://hicfiles.s3.amazonaws.com/hiseq/huvec/in-situ/combined.hic';...
%                 'https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic';...
%                 'http://hicfiles.s3.amazonaws.com/hiseq/rpe1/DarrowHuntley-2015/WT-combined.hic'};
    cellType = {'hHFFc6';...
                'hESC';...
                'IMR90';...
                'HUVEC';...
                'GM12878';...
                'RPE'};
    hicData = table(cellType,Data_Loc);

    % Hi-C matrix parameters
    if ~exist('chrSelect','var')||isempty(chrSelect);chrSelect = 14;end
    if ~exist('bpFrag','var')||isempty(bpFrag);bpFrag='BP';end
    if ~exist('binSize','var')||isempty(binSize);binSize=1E5;end

    H = cell(height(hicData));

    % load Hi-C data
    for iCt = 1:height(hicData)
        H{iCt} = hic2mat('oe','kr',hicData.Data_Loc{iCt},chrSelect,chrSelect,bpFrag,binSize,1,0);
    end

    %% Process data
    % remove disconnected nodes (genomic regions with a low number of aligned
    % reads)

    % [HTrim,badLocs] = hicTrim(cat(3,hHFFc6,hESC),2,.1);
    HTrim = cell(height(hicData));
    for iCt = 1:height(hicData)
        [HTrim{iCt},badLocs] = hicTrim(H{iCt},2,.1);
    end

    % compute the log2, change -inf values to minimum
    for iCt = 1:height(hicData)
        HTrim{iCt} = log2(HTrim{iCt});
        tempH = HTrim{iCt};
        tempH(tempH==-inf) = min(tempH(isfinite(tempH)));
        HTrim{iCt} = tempH;
    end

    %% VNE Computation
    % set VNE parameters
    preProcess = 'corr';

    % calculate VNE
    vnEntropy = zeros(height(hicData),1);
    for iCt = 1:height(hicData)
        vnEntropy(iCt) = hicVnEntropy(HTrim{iCt},[],[],preProcess);
    end

    %% Visualize Matrices with VNE
    % figure

    hicCMap = 1-((1-redblue(100))*.7);
    for iCt = 1:height(hicData)
        figure('name',hicData.cellType{iCt},'position',[50 50 700 500])

        imagesc(HTrim{iCt}), axis square
        title(sprintf('%s, VNE: %.2f', hicData.cellType{iCt}, vnEntropy(iCt)))
        colormap(hicCMap), caxis([-2 2]), colorbar
        ylabel('log_2(O/E)')

        set(get(gcf,'children'),'linewidth',2,'fontsize',20)
    end

    %% Save figures
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
      FigHandle = FigList(iFig);
      FigName   = get(FigHandle, 'Name');
      savefig(FigHandle, [Folder_Result, '\',FigName, '.fig']);
    end
end

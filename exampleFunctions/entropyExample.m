function [] = entropyExample(Data_Loc, Folder_Result, chrSelect, bpFrag, binSize)
    % Quantifying Structural Order through Network Entropy 
    % This example shows how entropy can be used to quantify
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
    %   * entropyExample.m created

    %% Set default Entropy parameters
    if ~exist('chrSelect','var')||isempty(chrSelect);chrSelect = 14;end
    if ~exist('bpFrag','var')||isempty(bpFrag);bpFrag='BP';end
    if ~exist('binSize','var')||isempty(binSize);binSize=1E5;end

    % load Hi-C data
    hHFFc6 = hic2mat('oe','kr',Data_Loc{1},chrSelect,chrSelect,bpFrag,binSize,1,0);
    hESC = hic2mat('oe','kr',Data_Loc{2},chrSelect,chrSelect,bpFrag,binSize,1,0);

    %% Process data
    % remove disconnected nodes (genomic regions with a low number of aligned
    % reads)
    [HTrim,badLocs] = hicTrim(cat(3,hHFFc6,hESC),2,.1);

    % compute the log2, change -inf values to minimum
%     HTrim = log2(HTrim);
%     for iA = 1:length(Data_Loc)
%         tempH = HTrim(:,:,iA);
%         tempH(tempH==-inf) = min(tempH(isfinite(tempH)));
% %         tempH = corr(tempH);
% %         tempH(isnan(tempH))=0;
%         HTrim(:,:,iA) = tempH;
%     end

    %% Entropy Computation
    % set Entropy parameters
    preProcess = 'corr';

    % calculate Entropy
    entropy = hicEntropy(HTrim,[],[],preProcess);

    %% Visualize Matrices with Entropy
    figure('name','Entropy Example', 'position',[50 50 1300 500])
    hicCMap = 1-((1-redblue(100))*.7);

    subplot(1,2,1)
    imagesc(HTrim(:,:,1)), axis square
    title(sprintf('Fibroblast, Entropy: %.2f',entropy(1)))
    colormap(hicCMap), caxis([-2 2]), colorbar
    ylabel('log_2(O/E)')

    subplot(1,2,2)
    imagesc(HTrim(:,:,2)), axis square
    title(sprintf('Embryonic Stem Cell, Entropy: %.2f',entropy(2)))
    colormap(hicCMap), caxis([-2 2]), colorbar

    set(get(gcf,'children'),'linewidth',2,'fontsize',20)
    linkaxes(get(gcf,'children'))

    %% Visualize Matrices with Entropy - corr
    % figure
%     figure('position',[50 50 1300 500])
%     hicCMap = 1-((1-redblue(100))*.7);
% 
%     subplot(1,2,2)
%     imagesc(corr(HTrim(:,:,1))), axis square
%     caxis([-1 1])
% 
%     subplot(1,2,1)
%     imagesc(corr(HTrim(:,:,2))), axis square
%     caxis([-1 1]), colorbar
%     ylabel('corr(log_2(O/E))')
% 
%     set(get(gcf,'children'),'linewidth',2,'fontsize',20)
%     linkaxes(get(gcf,'children'))
    
    %% Save figures
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
      FigHandle = FigList(iFig);
      FigName   = get(FigHandle, 'Name');
      savefig(FigHandle, [Folder_Result, '\',FigName, '.fig']);
      saveas(FigHandle, [Folder_Result, '\',FigName, '.png']);
    end
    
end
function [] = entropyExampleExpanded(Data_Loc, Folder_Result, chrSelect, bpFrag, binSize)
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

    HTrim = cell(height(hicData),1);
    for iCt = 1:height(hicData)
        [HTrim{iCt},badLocs] = hicTrim(H{iCt},2,.1);
    end
    
    
   for iCt = 1:height(hicData)
        eigvals{iCt} = eig(HTrim{iCt});
        eigvals{iCt} = sort(eigvals{iCt},'descend');
    end
    
    figure
    for iCt = 1:height(hicData)
        plot(eigvals{iCt})
        hold on
    end
    legend(cellType)

    % compute the log2, change -inf values to minimum
    for iCt = 1:height(hicData)
%         HTrim{iCt} = log2(HTrim{iCt});
%         tempH = HTrim{iCt};
%         tempH(tempH==-inf) = min(tempH(isfinite(tempH)));
        tempH = (tempH >= .95);
        tempH(tempH < .95) = 0;
        HTrim{iCt} = tempH;
    end

    %% VNE Computation
    % set VNE parameters
    preProcess = 'laplacian';

    % calculate VNE
    vnEntropy = zeros(height(hicData),1);
    vnEntropy2 = zeros(height(hicData),1);
    vnEntropyApprox = zeros(height(hicData),1);

    for iCt = 1:height(hicData)
        vnEntropy(iCt) = hicVnEntropy(HTrim{iCt},[],[],preProcess);
        vnEntropy2(iCt) = graph_entropy_test(HTrim{iCt},1);
        vnEntropyApprox(iCt) = hicVnEntropyApprox(HTrim{iCt});
    end

    %% Visualize Matrices with VNE
    hicCMap = 1-((1-redblue(100))*.7);
    for iCt = 1:height(hicData)
        figure('name',hicData.cellType{iCt},'position',[50 50 700 500])

        imagesc(HTrim{iCt}), axis square
        title(sprintf('%s, VNE: %.2f', hicData.cellType{iCt}, vnEntropy(iCt)))
        colormap(hicCMap), caxis([-2 2]), colorbar
        ylabel('log_2(O/E)')

        set(get(gcf,'children'),'linewidth',2,'fontsize',20)
    end
    
    %% Plot VNE over time  
    figure('Name','VNE and AVNE','position',[50 50 700 500])
    t_real = [0 2 5 7 15 80];
    plot(t_real,zscore(vnEntropy),'LineStyle','-','Color',rgb('darkgrey'),'Marker','o','MarkerEdgeColor','r', 'MarkerFaceColor','r','MarkerSize',6,'LineWidth',1.5); hold on;
    plot(t_real,zscore(vnEntropy2),'LineStyle','-','Color',rgb('darkgrey'),'Marker','o','MarkerEdgeColor','r', 'MarkerFaceColor','r','MarkerSize',6,'LineWidth',1.5); hold on;
    plot(t_real,zscore(vnEntropyApprox),'LineStyle','-','Color','k','Marker','^','MarkerEdgeColor','g', 'MarkerFaceColor','g','MarkerSize',6,'LineWidth',1.5); hold on;
    for t = 1:length(t_real)
        t_real_str{t} = sprintf('%d',t_real(t));
    end
    set(gca,'xtick',t_real) 
    set(gca,'xticklabel',t_real_str); 
    ylabel('Zscore'); 
    xlabel('Time (Day)');
    legend({'VNE', 'AVNE'});
    box off; legend boxoff ;

    %% Save figures
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
      FigHandle = FigList(iFig);
      FigName   = get(FigHandle, 'Name');
      savefig(FigHandle, [Folder_Result, '\',FigName, '.fig']);
      saveas(FigHandle, [Folder_Result, '\',FigName, '.png']);
    end
    
end

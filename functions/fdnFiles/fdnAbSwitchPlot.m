function [h] = fdnAbSwitchPlot(dataInfo,H,R,params)
%fdnAbSwitchPlot plots all sample A/B compartments with switch locations
%highlighted
%   
%   Inputs
%   dataInfo:   data structure with experiment sample info
%   H:          data structure containing all Hi-C data
%   R:          data structure containing all RNA-seq data
%   partParam:  data structure containing partitioning parameters
%   
%   Outputs
%   h:          Figure handles
%   
%   Version 1.0 (03/15/19)
%   Written by: Scott Ronquist
%   Contact: 	scotronq@umich.edu
%   Created: 	03/15/19
%   
%   Revision History:
%   v1.0 (03/15/19)
%   * fdnAbSwitchPlot.m created

%% set default parameters


%% plot
h = cell(length(params.numChr),1);
count = 1;
for iChr = params.numChr
    tempVar = var(H.s100kb.ABcomp{iChr},[],2);
    highVarLocs = tempVar > prctile(tempVar,params.varPerc);
    changeLocs = ~(all(H.s100kb.ABcomp{iChr}<0,2) | all(H.s100kb.ABcomp{iChr}>0,2));
    roiLocs = all([highVarLocs,changeLocs],2);
    
    % figure
    h{count} = figure;
    for iSamp = 1:size(H.s100kb.ABcomp{iChr},2)
        
        subplot(size(H.s100kb.ABcomp{iChr},2),1,iSamp)
        
        % add roiLocs background
        yLims = [min(H.s100kb.ABcomp{iChr}(:)) max(H.s100kb.ABcomp{iChr}(:))];
        
        bar(find(roiLocs),ones(sum(roiLocs),1)*yLims(1),'c',...
            'FaceColor','flat','EdgeColor','none'), hold on
        bar(find(roiLocs),ones(sum(roiLocs),1)*yLims(2),'c',...
            'FaceColor','flat','EdgeColor','none')
        
        % plot A/B bar
        b = bar(H.s100kb.ABcomp{iChr}(:,iSamp),'FaceColor','flat','EdgeColor','none');
        b.CData(H.s100kb.groupIdx{iChr}(:,iSamp)==1,:) = repmat([1 0 0],...
            sum(H.s100kb.groupIdx{iChr}(:,iSamp)==1),1);
        b.CData(H.s100kb.groupIdx{iChr}(:,iSamp)==2,:) = repmat([0 1 0],...
            sum(H.s100kb.groupIdx{iChr}(:,iSamp)==2),1);
        
        % add gene names
        if 1==params.geneNames
            text(find(roiLocs),H.s100kb.ABcomp{iChr}(roiLocs,iSamp),...
                R.s100kb.geneTrim{iChr}(roiLocs))
        end
        
        % format figure
        if iSamp == 1
            title(sprintf('Fiedler vector, Var percentile: %i, Chr %i',params.varPerc,iChr))
        elseif iSamp == size(H.s100kb.ABcomp{iChr},2)
            xlabel('Genomic Location, 100kb')
        end
        
        ylabel(dataInfo.sampleInfo.sample{iSamp})
        axis tight
        xticks(0:min([500,size(H.s100kb.ABcomp{iChr},1)]):size(H.s100kb.ABcomp{iChr},1))
    end
    
    linkaxes(get(gcf,'children'))
    set(get(gcf,'children'),'linewidth',2, 'fontsize', 15)
    
    count = count + 1;
end

end

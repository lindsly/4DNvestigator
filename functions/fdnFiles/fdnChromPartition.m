function [H,R] = fdnChromPartition(dataInfo,H,R,partParam)
%fdnChromPartition partitions the chromatin into nodal domains based on
%Hi-C contact map
%
%   Input
%   dataInfo: data structure with experiment sample info
%   H: data structure containing all Hi-C data
%   R: data structure containing all RNA-seq data
%   partParam: data structure containing partitioning parameters
%
%   Output
%   H: same as input, now with fields for chromatin partitioning
%   R: same as input, now with fields for chromatin partitioning
%
%
%   Version 1.1 (03/15/19)
%   Written by: Scott Ronquist
%   Contact:    scotronq@umich.edu
%   Created:    1/21/19
%   
%   Revision History:
%   v1.0 (1/21/19)
%   * fdnChromPartition.m created
%   v1.1 (03/15/19)
%   * formatted preamble
%   * added chr centromere finder (centLocStartTrim)

%% set parameter defaults
if ~isfield(partParam,'method')||isempty(partParam.method)
    partParam.method = 'fiedler';
end
if ~isfield('partParam','rnaSeqNorm')||isempty(partParam.rnaSeqNorm)
    partParam.rnaSeqNorm = [];
end
if ~isfield('partParam','chrDivide')||isempty(partParam.chrDivide)
    partParam.chrDivide = 'no';
end
if ~isfield('partParam','plotFlag')||isempty(partParam.plotFlag)
    partParam.plotFlag = 0;
end
if ~isfield('partParam','fdvSplitMethod')||isempty(partParam.fdvSplitMethod)
    partParam.fdvSplitMethod = 'sign';
end

%% get chr information from hic header
chrInfo = dataInfo.hicHeader.Chromosomes;
numChr = height(chrInfo);
numHicSamps = length(find(ismember(dataInfo.sampleInfo.dataType,'hic')));

%% Chromatin Partitioning
H.s100kb.ABcomp = cell(numChr,1);
H.s100kb.groupIdx = cell(numChr,1);
H.s100kb.centLocStartTrim = zeros(numChr,2);                                % matrix with start of chr centromere, and flag for whether chr was split A/B compartments
for iChr = 1:numChr
    % estimate centromere location
    idx = 1:size(H.s100kb.oe{iChr},1);
    idxTrim = idx(~H.s100kb.oeTrimBadLocs{iChr});
    [temp,H.s100kb.centLocStartTrim(iChr,1)] = max(diff(idxTrim));
    
    % decide whether to use chr split
    % if:   centromere region is large, and near chr middle
    % then: calculate A/B comp on split chr
    if temp > length(idx)*.01 && H.s100kb.centLocStartTrim(iChr)/length(idx) < .8 &&...
             H.s100kb.centLocStartTrim(iChr)/length(idx) > .2
        partParam.chrDivide = 'yes';
        H.s100kb.centLocStartTrim(iChr,2) = 1;
    else
        partParam.chrDivide = 'no';
        H.s100kb.centLocStartTrim(iChr,2) = 0;
    end
    
    % get Chromatin Partitioning
    for iSample = 1:numHicSamps
        fprintf('calculating A/B compartments, Sample: (%d/%d), chr:%d...\n',iSample,numHicSamps,iChr)
        if ~isempty(H.s100kb.oeTrim{iChr}(:,:,iSample))
            
            hTemp = H.s100kb.oeTrim{iChr}(:,:,iSample);
            hTemp(hTemp>prctile(hTemp(:),99)) = prctile(hTemp(:),99);
            
            rTemp = log2(R.s100kb.tpmMeanTrim{iChr}(:,iSample)+1);
            [H.s100kb.ABcomp{iChr}(:,iSample),H.s100kb.groupIdx{iChr}(:,iSample)] =...
                hicABcomp(hTemp,partParam.method,rTemp,partParam.rnaSeqNorm,...
                partParam.chrDivide,partParam.plotFlag,partParam.fdvSplitMethod,...
                H.s100kb.centLocStartTrim(iChr,1));
        end
    end
end

% %% plot each partition
% if partParam.plotFlag
%     fn = [dataInfo.path.output,dataInfo.delim,'figures',dataInfo.delim,'chromPart'];
%     mkdir(fn)
%     for iChr = 1:numChr
%         for iSample = 1:numHicSamps
%             fprintf('plotting figure, Sample: (%d/%d), chr:%d...\n',iSample,numHicSamps,iChr)
%             
%             % skip if doesn't exist
%             try
%                 temp = H.s100kb.ABcomp{iChr}(:,iSample);
%             catch
%                 continue
%             end
%             
%             % get figure properties
%             hicCMap = [ones(64,1),[1:-1/63:0]',[1:-1/63:0]'];
%             figure('position',[100 100 700 1000])
%             
%             % plot RNA-seq
%             ax1 = subplot(6,1,1);
%             rTemp = log2(R.s100kb.tpmMeanTrim{iChr}(:,iSample)+1);
%             bar(rTemp), axis tight
%             title(sprintf('%s, T:%i, Chr%s',dataInfo.sampleInfo.sample{iSample},...
%                 dataInfo.sampleInfo.timePoint(iSample),chrInfo.chr{iChr}))
%             ylabel('log_2 TPM')
%             
%             % plot ABcomp
%             ax2 = subplot(6,1,2);
%             tempAB = H.s100kb.ABcomp{iChr}(:,iSample);
%             b = bar(tempAB,'FaceColor','flat','EdgeColor','none');
%             b.CData(H.s100kb.groupIdx{iChr}(:,iSample)==1,:) = repmat([1 0 0],...
%                 sum(H.s100kb.groupIdx{iChr}(:,iSample)==1),1);
%             b.CData(H.s100kb.groupIdx{iChr}(:,iSample)==2,:) = repmat([0 1 0],...
%                 sum(H.s100kb.groupIdx{iChr}(:,iSample)==2),1);
%             axis tight
%             ylabel(partParam.method)
%             
%             % plot Hi-C
%             ax3 = subplot(6,1,3:6);
%             hTemp = H.s100kb.oeTrim{iChr}(:,:,iSample);
%             climMain = [0 prctile(hTemp(:),90)];
%             hTemp(hTemp>prctile(hTemp(:),99)) = prctile(hTemp(:),99);
%             imagesc(hTemp), axis square
%             colormap(ax3,hicCMap); %colorbar
%             caxis(climMain)
%             ylabel(sprintf('chr%i Hi-C map',iChr))
%             
%             % figure format
%             set(get(gcf,'children'),'linewidth',2,'fontsize',15)
%             linkaxes(get(gcf,'children'),'x')
%             
%             % save figure
%             saveas(gcf,sprintf('%s%ss%s_t%i_chr%s.fig',fn,...
%                 dataInfo.delim,dataInfo.sampleInfo.sample{iSample},...
%                 dataInfo.sampleInfo.timePoint(iSample),...
%                 chrInfo.chr{iChr}))
%             
%             close all
%         end
%     end
% end
% 
% %% add average A/B to gene table
% if 1==0
%     % create gene AB table
%     R.abTable = R.TPM(:,1:6);
%     ab = zeros(height(R.abTable),size(H.s100kb.ABcomp{1},2));
%     
%     % loop through gene list and find AB for each gene
%     for iGene = 1:height(R.abTable)
%         disp(iGene/height(R.abTable))
%         
%         % find gene name in trim 100kb bin genenames
%         temp = mean(H.s100kb.ABcomp{R.abTable.chr(iGene)}...
%             (cellfun(@(x) any(strcmp(x,R.abTable.geneName{iGene})),...
%             R.s100kb.geneTrim{R.abTable.chr(iGene)}),:),1);
%         
%         ab(iGene,:) = temp;
%     end
%     
%     R.abTable = [R.abTable,table(ab)];
% end
%% find AB switch region



end

function [H] = gsfatChromPartition(dataInfo,H,R)
%gsfatChromPartition partitions the chromatin into nodal domains
%   Detailed explanation goes here

numChr = height(dataInfo.chrSizes);
numHicSamps = length(find(ismember(dataInfo.sampleInfo.dataType,'hic')));

%% Chromatin Partitioning
%%% PARAMETERS vvv
hicParam.method='fiedler';
hicParam.rnaSeq=[];
hicParam.rnaSeqNorm=[];
hicParam.chrDivide='no';
hicParam.plotFlag=0;
hicParam.fdvSplitMethod='sign';
%%% PARAMETERS ^^^

H.s100kb.ABcomp = cell(numChr,1);
H.s100kb.groupIdx = cell(numChr,1);
for iChr = 1:numChr
    for iSample = 1:numHicSamps
        fprintf('AB comp, Sample: (%d/%d), chr:%d...\n',iSample,numHicSamps,iChr)
        
        % get Chromatin Partitioning
        hTemp = H.s100kb.oeTrim{iChr}(:,:,iSample);
        [H.s100kb.ABcomp{iChr}(:,iSample),H.s100kb.groupIdx{iChr}(:,iSample)] =...
            hicABcomp(hTemp,hicParam.method,hicParam.rnaSeq,hicParam.rnaSeqNorm,...
            hicParam.chrDivide,hicParam.plotFlag,hicParam.fdvSplitMethod);
        
        % fix Chromatin Partitioning sign by correlating with # of genes
        tempGeneNum = cellfun(@length,R.s100kb.geneTrim{iChr});
        if mean(tempGeneNum(H.s100kb.groupIdx{iChr}(:,iSample)==1)) >...
                mean(tempGeneNum(H.s100kb.groupIdx{iChr}(:,iSample)==2))
            H.s100kb.ABcomp{iChr}(:,iSample) = -H.s100kb.ABcomp{iChr}(:,iSample);
            H.s100kb.groupIdx{iChr}(:,iSample) = -H.s100kb.groupIdx{iChr}(:,iSample)+3;
        end
    end
end


end


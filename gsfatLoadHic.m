function [H] = gsfatLoadHic(dataInfo)
%gsfatLoadHiC Loads and formats Hi-C data specified in dataInfo
%   Detailed explanation goes here

numChr = height(dataInfo.chrSizes);

%% get Hi-C sample Locations
sampleFn = dataInfo.sampleInfo.path(ismember(dataInfo.sampleInfo.dataType,'hic'));

%% get chr sizes to stitch together
chrSizesSorted = zeros(numChr,1);
for iChr = 1:numChr
    if iChr == numChr+1
        iChr_ = 'X';
    else
        iChr_ = num2str(iChr);
    end
    for i = 1:height(dataInfo.chrSizes)
        if strcmp(dataInfo.chrSizes{i,1}{1}(4:end),iChr_)
            chrSizesSorted(iChr) = dataInfo.chrSizes{i,2};
        end
    end
end

%% get hic header info
[~,dataInfo.hicHeader] = hic2mat('observed','NONE',sampleFn{1},...
    numChr,numChr,'BP',1E6,hicParam.intraFlag);

%% load Hi-C data - 100kb
%%% PARAMETERS vvv
hicParam.binType = 'BP';
hicParam.binSize = 1E5;
hicParam.norm1d = 'KR';
hicParam.norm3d = 'oe';
hicParam.intraFlag = 1;
%%% PARAMETERS ^^^

chrBinSizes = ceil(dataInfo.chrSizes{:,2}/hicParam.binSize);
H.s100kb.oe = cell(numChr,1);
H.s100kb.kr = cell(numChr,1);

for iChr = 1:numChr
    H.s100kb.oe{iChr} = zeros(chrBinSizes(iChr),chrBinSizes(iChr),length(sampleFn));
    H.s100kb.kr{iChr} = zeros(chrBinSizes(iChr),chrBinSizes(iChr),length(sampleFn));
    for iSample = 1:length(sampleFn)
        fprintf('loading 100kb Hi-C. Sample: (%d/%d), chr:%d...\n',iSample,length(sampleFn),iChr)
        
        temp = hic2mat(hicParam.norm3d,hicParam.norm1d,sampleFn{iSample},...
            iChr,iChr,hicParam.binType,hicParam.binSize,hicParam.intraFlag);
        temp(isnan(temp)) = 0;
        temp = max(cat(3,temp,temp'),[],3);
        H.s100kb.oe{iChr}(1:length(temp),1:length(temp),iSample) = temp;
        
        temp = hic2mat('observed',hicParam.norm1d,sampleFn{iSample},...
            iChr,iChr,hicParam.binType,hicParam.binSize,hicParam.intraFlag);
        temp(isnan(temp)) = 0;
        temp = max(cat(3,temp,temp'),[],3);
        H.s100kb.kr{iChr}(1:length(temp),1:length(temp),iSample) = temp;
    end
end

%% load Hi-C data - 1Mb
%%% PARAMETERS vvv
hicParam.binSize = 1E6;
%%% PARAMETERS ^^^

H.s1mb.chrStart = [1;cumsum(ceil(chrSizesSorted/hicParam.binSize))+1];

% get Hi-C data
H.s1mb.oe = zeros(H.s1mb.chrStart(end)-1,H.s1mb.chrStart(end)-1,length(sampleFn));
H.s1mb.kr = zeros(H.s1mb.chrStart(end)-1,H.s1mb.chrStart(end)-1,length(sampleFn));
for iSample = 1:length(sampleFn)
    for iChr1 = 1:numChr
        for iChr2 = iChr1:numChr
            iChr1_ = iChr1;iChr2_ = iChr2;
            if iChr1 == numChr+1;  iChr1_ = 'X';end
            if iChr2 == numChr+1;  iChr2_ = 'X';end
            fprintf('loading 1Mb Hi-C. Sample: (%d/%d), chr1:%d, chr2:%d...\n',...
                iSample,length(sampleFn),iChr1,iChr2)
            
            temp = hic2mat(hicParam.norm3d,hicParam.norm1d,sampleFn{iSample},...
                iChr1_,iChr2_,hicParam.binType,hicParam.binSize);
            
            if chrSizesSorted(iChr2)>chrSizesSorted(iChr1); temp=temp';end
            H.s1mb.oe(H.s1mb.chrStart(iChr1):H.s1mb.chrStart(iChr1)+size(temp,1)-1,...
                H.s1mb.chrStart(iChr2):H.s1mb.chrStart(iChr2)+size(temp,2)-1,iSample) = temp;
            
            temp = hic2mat('observed',hicParam.norm1d,sampleFn{iSample},...
                iChr1_,iChr2_,hicParam.binType,hicParam.binSize);
            
            if chrSizesSorted(iChr2)>chrSizesSorted(iChr1); temp=temp';end
            H.s1mb.kr(H.s1mb.chrStart(iChr1):H.s1mb.chrStart(iChr1)+size(temp,1)-1,...
                H.s1mb.chrStart(iChr2):H.s1mb.chrStart(iChr2)+size(temp,2)-1,iSample) = temp;
        end
    end
    H.s1mb.oe(:,:,iSample) = max(cat(3,H.s1mb.oe(:,:,iSample),H.s1mb.oe(:,:,iSample)'),[],3);
    H.s1mb.kr(:,:,iSample) = max(cat(3,H.s1mb.kr(:,:,iSample),H.s1mb.kr(:,:,iSample)'),[],3);
end

%% trim data (remove unmappable regions) - 100kb and 1Mb
%%% PARAMETERS vvv
hicParam.numDiag = 2;
hicParam.numSparse = .1;
%%% PARAMETERS ^^^

% trim 100kb all together
H.s100kb.oeTrim = cell(numChr,1);
H.s100kb.krTrim = cell(numChr,1);
H.s100kb.genesTrim = cell(numChr,1);
H.s100kb.oeTrimBadLocs = cell(numChr,1);
for iChr = 1:numChr
    fprintf('Trimming Hi-C Chr:%d...\n',iChr)
    [H.s100kb.oeTrim{iChr},H.s100kb.oeTrimBadLocs{iChr}] = hicTrim(H.s100kb.oe{iChr},...
        hicParam.numDiag,hicParam.numSparse);
    H.s100kb.krTrim{iChr} = H.s100kb.kr{iChr}(~H.s100kb.oeTrimBadLocs{iChr},~H.s100kb.oeTrimBadLocs{iChr},:);
end

% trim 1Mb all together
[H.s1mb.oeTrim,H.s1mb.oeTrimBadLocs] = hicTrim(H.s1mb.oe,...
    hicParam.numDiag,hicParam.numSparse);
H.s1mb.krTrim = H.s1mb.kr(~H.s1mb.oeTrimBadLocs,~H.s1mb.oeTrimBadLocs,:);

% update chr Start site with trimmed Locs
H.s1mb.chrStartTrim = H.s1mb.chrStart;
for iChr = 1:numChr
    H.s1mb.chrStartTrim(iChr) = H.s1mb.chrStart(iChr)-...
        sum(H.s1mb.oeTrimBadLocs(1:H.s1mb.chrStart(iChr)-1));
end

end


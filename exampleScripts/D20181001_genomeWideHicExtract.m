% goal: write a script to extract genome

%%% PARAMETERS vvv
hicParam.binSize = 1E6;
hicParam.binType = 'BP';
hicParam.norm1d = 'KR';
hicParam.norm3d = 'oe';
hicParam.intraFlag = 1;
%%% PARAMETERS ^^^

% get chr sizes to stitch together
chrSizes = readtable('hg19.chrom.sizes','filetype','text');

H.s1mb.chrStart = [1;cumsum(ceil(chrSizes{:,1:23}./hicParam.binSize))+1];

% get Hi-C data
H.s1mb.oe = zeros(H.s1mb.chrStart(end)-1,H.s1mb.chrStart(end)-1,length(samplesNumber));
for iSample = 1:length(samplesNumber)
    for iChr1 = 1:23
        for iChr2 = iChr1:23
            iChr1_ = iChr1;iChr2_ = iChr2;
            if iChr1 == 23;  iChr1_ = 'X';end
            if iChr2 == 23;  iChr2_ = 'X';end
            fprintf('loading 1Mb Hi-C. Sample: (%d/%d), chr1:%d, chr2:%d...\n',...
                iSample,length(samplesNumber),iChr1,iChr2)
            
            temp = hic2mat(hicParam.norm3d,hicParam.norm1d,...
                [dataFolder,'\Sample_',num2str(samplesNumber(iSample)),'\inter_30.hic'],...
                iChr1_,iChr2_,hicParam.binType,hicParam.binSize);
            
        end
    end
    H.s1mb.oe(:,:,iSample) = max(cat(3,H.s1mb.oe(:,:,iSample),H.s1mb.oe(:,:,iSample)'),[],3);
end
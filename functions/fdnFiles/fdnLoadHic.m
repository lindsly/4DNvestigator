function [H] = fdnLoadHic(dataInfo)
%fdnLoadHic Loads and formats Hi-C data specified in dataInfo
%
%   Input
%   dataInfo:   structure containing all experiment metadata
%
%   Output
%   H:          structure containing all Hi-C data needed for 4DNvestigator
%
%   Version 1.1 (03/15/19)
%   Written by: Scott Ronquist
%   Contact:    scotronq@umich.edu
%   Created:    12/13/18
%   
%   Revision History:
%   v1.0 (12/13/18)
%   * fdnLoadHic.m created
%   v1.1 (03/15/19)
%   * formatted preamble

%% get information dataInfo
% get chr information from hic header
chrInfo = dataInfo.hicHeader.Chromosomes;
numChr = height(chrInfo);

% get Hi-C sample Locations
sampleFn = dataInfo.sampleInfo.path(ismember(dataInfo.sampleInfo.dataType,'hic'));

% Create the waitbar and determine intialization properties
waitBar = waitbar(0,'loading 100kb Hi-C...','Name','Loading Data');
set(findall(waitBar),'Units', 'normalized');    % Change the Units Property of the figure and all the children
set(waitBar,'Position', [0.25 0.4 0.5 0.08]);   % Change the size of the figure
totalWait = numChr*length(sampleFn);
currentWait = 1;
avTime = 0;

%% load Hi-C data - 100kb
% this loads intra-chr Hi-C
%%% PARAMETERS vvv
hicParam.binType = 'BP';
hicParam.binSize = 1E5;
hicParam.norm1d = 'KR';
hicParam.norm3d = 'oe';
hicParam.intraFlag = 1;
%%% PARAMETERS ^^^

chrBinSizes = ceil(chrInfo.chrLength/hicParam.binSize);
H.s100kb.oe = cell(numChr,1);
H.s100kb.kr = cell(numChr,1);

for iChr = 1:numChr
    H.s100kb.oe{iChr} = zeros(chrBinSizes(iChr),chrBinSizes(iChr),length(sampleFn));
    H.s100kb.kr{iChr} = zeros(chrBinSizes(iChr),chrBinSizes(iChr),length(sampleFn));
    for iSample = 1:length(sampleFn)
        % for estimation of time
        if iChr==1 && iSample==1;tic;end
        
        % update load bar
        waitbar(currentWait/totalWait,waitBar,...
            sprintf('Loading 100kb Hi-C. Sample: (%d/%d), chr:%s, estimated time: %i secs...',...
            iSample,length(sampleFn),chrInfo.chr{iChr},round(avTime*(totalWait-currentWait))));
        
        % extract o/e normalized
        temp = hic2mat(hicParam.norm3d,hicParam.norm1d,sampleFn{iSample},...
            chrInfo.chr{iChr},chrInfo.chr{iChr},hicParam.binType,hicParam.binSize,hicParam.intraFlag);
        temp(isnan(temp)) = 0;
        temp = max(cat(3,temp,temp'),[],3);
        H.s100kb.oe{iChr}(1:length(temp),1:length(temp),iSample) = temp;
        
        % extract observed as well
        temp = hic2mat('observed',hicParam.norm1d,sampleFn{iSample},...
            chrInfo.chr{iChr},chrInfo.chr{iChr},hicParam.binType,hicParam.binSize,hicParam.intraFlag);
        temp(isnan(temp)) = 0;
        temp = max(cat(3,temp,temp'),[],3);
        H.s100kb.kr{iChr}(1:length(temp),1:length(temp),iSample) = temp;
        
        % update load bar variable
        currentWait = currentWait+1;
        if iChr==1 && iSample==1;avTime = toc;end % for estimation of time to load
    end
    
end
close(waitBar)

%% load Hi-C data - 1Mb
% this loads both intra- and inter-chr Hi-C
%%% PARAMETERS vvv
hicParam.binSize = 1E6;
%%% PARAMETERS ^^^

% Create the waitbar and determine intialization properties
waitBar = waitbar(0,'loading 1Mb Hi-C...','Name','Loading Data');
set(findall(waitBar),'Units', 'normalized');    % Change the Units Property of the figure and all the children
set(waitBar,'Position', [0.25 0.4 0.5 0.08]);   % Change the size of the figure
totalWait = (nchoosek(numChr,2)+numChr)*length(sampleFn);
currentWait = 1;
avTime = 0;

H.s1mb.chrStart = [1;cumsum(ceil(chrInfo.chrLength/hicParam.binSize))+1];

% get Hi-C data
H.s1mb.oe = zeros(H.s1mb.chrStart(end)-1,H.s1mb.chrStart(end)-1,length(sampleFn));
H.s1mb.kr = zeros(H.s1mb.chrStart(end)-1,H.s1mb.chrStart(end)-1,length(sampleFn));
for iSample = 1:length(sampleFn)
    for iChr1 = 1:numChr
        for iChr2 = iChr1:numChr
            % for estimation of time
            if iChr1==1 && iChr2==1 && iSample==1;tic;end
            
            % update load bar
            waitbar(currentWait/totalWait,waitBar,...
                sprintf('Loading 1Mb Hi-C. Sample: (%d/%d), chr1:%s, chr2:%s, estimated time: %i secs...',...
                iSample,length(sampleFn),chrInfo.chr{iChr1},chrInfo.chr{iChr2},...
                round(avTime*(totalWait-currentWait))));
            
            % o/e
            temp = hic2mat(hicParam.norm3d,hicParam.norm1d,sampleFn{iSample},...
                chrInfo.chr{iChr1},chrInfo.chr{iChr2},hicParam.binType,hicParam.binSize);
            %if chrSizesSorted(iChr2)>chrSizesSorted(iChr1); temp=temp';end
            H.s1mb.oe(H.s1mb.chrStart(iChr1):H.s1mb.chrStart(iChr1)+size(temp,1)-1,...
                H.s1mb.chrStart(iChr2):H.s1mb.chrStart(iChr2)+size(temp,2)-1,iSample) = temp;
            
            % observed
            temp = hic2mat('observed',hicParam.norm1d,sampleFn{iSample},...
                chrInfo.chr{iChr1},chrInfo.chr{iChr2},hicParam.binType,hicParam.binSize);
            %if chrSizesSorted(iChr2)>chrSizesSorted(iChr1); temp=temp';end
            H.s1mb.kr(H.s1mb.chrStart(iChr1):H.s1mb.chrStart(iChr1)+size(temp,1)-1,...
                H.s1mb.chrStart(iChr2):H.s1mb.chrStart(iChr2)+size(temp,2)-1,iSample) = temp;
            
            % update load bar variable
            currentWait = currentWait+1;
            
            % for estimation of time to load
            if iChr1==1 && iChr2==1 && iSample==1;avTime = toc;end
        end
    end
    H.s1mb.oe(:,:,iSample) = max(cat(3,H.s1mb.oe(:,:,iSample),H.s1mb.oe(:,:,iSample)'),[],3);
    H.s1mb.kr(:,:,iSample) = max(cat(3,H.s1mb.kr(:,:,iSample),H.s1mb.kr(:,:,iSample)'),[],3);
end
close(waitBar)

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
    
    % estimate centromere location
    idx = 1:size(H.s100kb.oe{iChr},1);
    idxTrim = idx(~H.s100kb.oeTrimBadLocs{iChr});
    
    [~,H.s100kb.centLocStart{iChr}] = max(diff(idxTrim));
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

%% get centromere loc

end


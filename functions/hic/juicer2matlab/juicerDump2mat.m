function [Hj] = juicerDump2mat(juicerFn,intraFlag,specificBins,dataType)
%juicer2mat creates a Hi-C adjacency matrix from "juicer dump" output
%   This function outputs a square adjacency matrix from juicer output.
%   juicer dump command creates a 3 column tab-delimited file: bin1 bin2
%   count
%
%   Input
%   juicerFn:      Juicer file name or juicer_tools dump output (N x 3 matrix)
%   intraFlag:      Denote intra-chr (1) vs inter-chr (0) (default: 0)
%   specificBins:  Specify bins to output (default: all bins output)
%   dataType:      Variable number type (default: double)
%
%   Output
%   Hj:             MATLAB formatted Hi-C matrix
%
%   Scott Ronquist, scotronq@umich.edu. 1/22/19

%% set default parameters
if ~exist('intraFlag','var')||isempty(intraFlag); intraFlag = 0; end
if ~exist('specificBins','var')||isempty(specificBins); specificBins = []; end
if ~exist('dataType','var')||isempty(dataType); dataType = 'double'; end

%% extract juicer txt data
if ischar(juicerFn)
    fileID = fopen(juicerFn);
    C = textscan(fileID,'%d %d %d');
    fclose(fileID);
    C = cell2mat(C);
    C(:,1:2) = C(:,1:2)+1;
else
    C = juicerFn;
    C(:,1:2) = C(:,1:2)+1;
end

%% extra
if ~isempty(specificBins)
    [~,~,bin1] = histcounts(C(:,1),[specificBins(:,1);specificBins(end,2)]);
    [~,~,bin2] = histcounts(C(:,2),[specificBins(:,1);specificBins(end,2)]);
    keep_idx = diff([bin1,bin2],1,2)==0;
    C(~keep_idx,:)=[];
end

%% create matrix
if intraFlag
    HSize = max(max(C(:,1:2)));
    
    Hj = zeros(HSize,dataType);%,'int32');
    idx = sub2ind(size(Hj),C(:,1),C(:,2));
    Hj(idx) = C(:,3);
    idx = sub2ind(size(Hj),C(:,2),C(:,1));
    Hj(idx) = C(:,3);
else
    HSize = max(C(:,1:2));
    
    Hj = zeros(HSize,dataType);%,'int32');
    idx = sub2ind(size(Hj),C(:,1),C(:,2));
    Hj(idx) = C(:,3);
end

% remove NaNs and Infs
Hj(isnan(Hj)) = 0;
Hj(isinf(Hj)) = 0;

end


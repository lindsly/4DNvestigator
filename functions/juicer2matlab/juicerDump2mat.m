function [H_j] = juicerDump2mat(juicerFn,intraFlag,specificBins,dataType)
%juicer2mat creates a Hi-C adjacency matrix from "juicer dump" output
%   This function outputs a square adjacency matrix from juicer output.
%   juicer dump command creates a 3 column tab-delimited file: bin1 bin2
%   count
%
%   juicer_fn: juicer file name or juicer_tools dump output (N x 3 matrix)
%   intraFlag: denote intra-chr (1) vs inter-chr (0) (default: 0)
%   specific_bins: specify bins to output (default: all bins output)
%   data_type: variable number type (default: double)
%
%   Scott Ronquist, 6/27/18

%% to be used later for whole genome input
if ~exist('intraFlag','var')||isempty(intraFlag); intraFlag = 0; end
if ~exist('specificBins','var')||isempty(specificBins); specificBins = []; end
if ~exist('dataType','var')||isempty(dataType); dataType = 'double'; end

% extract juicer txt data
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
% cheat fix for inter vs intra, probably should make more robust
if intraFlag
    H_size = max(max(C(:,1:2)));
    
    H_j = zeros(H_size,dataType);%,'int32');
    idx = sub2ind(size(H_j),C(:,1),C(:,2));
    H_j(idx) = C(:,3);
    idx = sub2ind(size(H_j),C(:,2),C(:,1));
    H_j(idx) = C(:,3);
else
    H_size = max(C(:,1:2));
    
    H_j = zeros(H_size,dataType);%,'int32');
    idx = sub2ind(size(H_j),C(:,1),C(:,2));
    H_j(idx) = C(:,3);
end

%% extra
% transpose is slower somehow
% tic
% H_j = zeros(H_size);
% idx = sub2ind(size(H_j),C(:,1),C(:,2));
% H_j(idx) = C(:,3);
% H_j = H_j+triu(H_j,1)';
% toc

end


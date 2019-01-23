function [rnaBin,rnaGeneNamesBin]= rna2bin(rnaseq,geneName,geneLocs,binSize,chrBinLengths,method)
%rna2bin converts RNA-seq data from gene to bin-level (ie 100kb, 1MB)
%
%   Input
%   rnaseq:         Rna-seq data, N x 1 vector
%   geneName:       Gene names, N x 1 vector
%   geneLocs:       Location of gene, N x 3 vector [chr start end]
%   binSize:        Bin resolution of output (eg 1E5 for 100kb)
%   chrBinLengths:  Length of chr, 22 x 1 vector (needed for match with Hi-C)
%   method:         'sum' or 'mean'
%
%   Output
%   rnaBin:         RNA-seq values, binned by "binSize"
%   rnaGeneNamesBin:Gene names within each bin
%
%   Scott Ronquist, 1/22/18

%% set default parameters
if nargin <6; method='sum';end

%%
for chr=1:length(chrBinLengths)
    fprintf('binning chr%d \n',chr)
    chrIdx = geneLocs(:,1)==chr;
    rnaBin{chr,1} = zeros(chrBinLengths(chr),size(rnaseq,2));
    
    chrRnaseq = rnaseq(chrIdx,:);
    chrGeneName = geneName(chrIdx,:);
    chrGeneLocs = geneLocs(chrIdx,:);
    
    %% overlap faction between gene and bin
    overlapRatio = bboxOverlapRatio([[1:binSize:chrBinLengths(chr)*binSize]',...
        zeros(chrBinLengths(chr),1), ones(chrBinLengths(chr),1)*binSize,...
        ones(chrBinLengths(chr),1)],...
        [chrGeneLocs(:,2), zeros(length(chrGeneLocs),1),...
        diff(chrGeneLocs(:,2:3),1,2), ones(length(chrGeneLocs),1)],'min');
    
    maxLengthRef = max([diff(chrGeneLocs(:,2:3),1,2),ones(size(chrGeneLocs,1),1)*binSize],[],2);
    geneOverlapPercentages = (overlapRatio./repmat(maxLengthRef',size(overlapRatio,1),1))*binSize;
    
    for i = 1:chrBinLengths(chr)
        rnaGeneNamesBin{chr,1}{i,1} = chrGeneName(logical(geneOverlapPercentages(i,:)));
    end
    
    % take either the mean or sum
    switch method
        case 'sum'
            for tp = 1:size(rnaseq,2)
                rnaBin{chr,1}(:,tp) = geneOverlapPercentages*chrRnaseq(:,tp);
            end
        case 'mean'
            error('mean method is depreciated')
    end
end

end
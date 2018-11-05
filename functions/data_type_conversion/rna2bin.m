function [rna_bin,rna_gene_names_bin]= rna2bin(rnaseq,gene_name,gene_locs,bin_size,chr_bin_lengths,method)
%rna2bin converts RNA-seq data from gene to bin-level (ie 100kb, 1MB)
%   rnaseq: rnaseq data, N x 1 vector
%   gene_name: gene names, N x 1 vector
%   gene_locs: location of gene, N x 3 vector [chr start end]
%   bin_size: bin resolution of output (eg 1E5 for 100kb)
%   chr_bin_lengths: length of chr, 22 x 1 vector (needed for mathc with Hi-C)
%   method: sum or mean
%
%   Scott Ronquist, 4/21/18

if nargin <6; method='sum';end

for chr=1:length(chr_bin_lengths)
    fprintf('binning chr%d \n',chr)
    chr_idx = gene_locs(:,1)==chr;
    rna_bin{chr,1} = zeros(chr_bin_lengths(chr),size(rnaseq,2));
    
    chr_rnaseq = rnaseq(chr_idx,:);
    chr_gene_name = gene_name(chr_idx,:);
    chr_gene_locs = gene_locs(chr_idx,:);
    
    %% overlap faction between gene and bin
    overlapRatio = bboxOverlapRatio([[1:bin_size:chr_bin_lengths(chr)*bin_size]',...
        zeros(chr_bin_lengths(chr),1), ones(chr_bin_lengths(chr),1)*bin_size,...
        ones(chr_bin_lengths(chr),1)],...
        [chr_gene_locs(:,2), zeros(length(chr_gene_locs),1),...
        diff(chr_gene_locs(:,2:3),1,2), ones(length(chr_gene_locs),1)],'min');
    
    max_length_ref = max([diff(chr_gene_locs(:,2:3),1,2),ones(size(chr_gene_locs,1),1)*bin_size],[],2);
    gene_overlap_percentages = (overlapRatio./repmat(max_length_ref',size(overlapRatio,1),1))*bin_size;
    
    for i = 1:chr_bin_lengths(chr)
        rna_gene_names_bin{chr,1}{i,1} = chr_gene_name(logical(gene_overlap_percentages(i,:)));
    end
    
    switch method
        case 'sum'
            for tp = 1:size(rnaseq,2)
                rna_bin{chr,1}(:,tp) = gene_overlap_percentages*chr_rnaseq(:,tp);
            end
        case 'mean'
            error('NOT FISNISHED, IN PROGRESS')
    end
end

end
function [R_table] = rnaseq_txt2mat(fn,genome)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 1; fn='GSE40705_FFPE_gene.txt';end
if nargin < 2; genome='hg19';end

%% load biomart data
switch genome
    case {'hg19','GRCh37'}
        biomart = readtable('mart_export_ensembl_hg37_info.txt');
    case {'hg38','GRCh38'}
        biomart = readtable('mart_export_ensembl_hg37_info.txt');
end

%% 
rnaseq_txt = readtable(fn);
genes = unique(rnaseq_txt.Var3);

ismember(rnaseq_txt.Var3,biomart.HGNCSymbol)


R_table
end


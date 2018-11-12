% This script performs all default analysis methods for the Genome
% Structure Function Analysis Toolbox (GSFAT)
% https://github.com/scotronq/srMatlabFunctions

% Variables
% dataInfo: MATLAB Structure that stores information that explains the input data.
% Includes data dependencies and file locations

% R: MATLAB Structure that contains all RNA-seq data
% H: MATLAB Structure that contains all Hi-C data

% Example

% Scott Ronquist, 10/16/18

%% default start
clear
close all
restoredefaultpath
addpath(genpath('.'))

isTest = 0;

%% Input data Information
if 1==isTest
    load('testDataInfo.mat')
else
    [dataInfo] = gsfatLoadUserInput;
end

%GOI list? gene Info table?

%% load Hi-C
% loads and formats Hi-C data output from Juicer
if 1==isTest
    load('testDataInfo.mat','H')
else
    H = gsfatLoadHic(dataInfo);
end

%% load RNA-seq
% loads and formats RNA-seq data output from RSEM
% H variable is loaded first and passed into this function to inform on
% which regions were "trimmed" (removed due to sparse coverage)
if 1==0
    load('testDataInfo.mat','R')
else
    R = gsfatLoadRnaseq(dataInfo,H);
end

%% Analysis
%% Chromatin partitioning
[H] = gsfatChromPartition(dataInfo,H,R);

%% RNA-seq differential expression


%% 4DN feature analyzer


%% matrix comparisons










%% load RNA-seq
% loads and formats RNA-seq data output from RSEM
R = rsem2mat(dataInfo.path.rnaseq,dataInfo.refGenome);

R_fields = fields(R);

for i = 1:length(R_fields)
    R.(R_fields{i})(~cellfun(@isnumeric,R.(R_fields{i}).chr),:) = [];
    R.(R_fields{i}).chr = cell2mat(R.(R_fields{i}).chr);
end

%100kb bins
chr_bin_lengths = diff(H.s100kb.chr_info{:,2:3},1,2)+1;
[R.s100kb.tpm,R.s100kb.gene] = rna2bin(R.TPM{:,7:end},R.TPM.gene_name,...
    [R.TPM.chr R.TPM.gene_start R.TPM.gene_end],1E5,chr_bin_lengths(1:22));

%1mb bins
chr_bin_lengths = diff(H.s1mb.chr_info{:,2:3},1,2)+1;
[R.s1mb.tpm,R.s1mb.gene] = rna2bin(R.TPM{:,7:end},R.TPM.gene_name,...
    [R.TPM.chr R.TPM.gene_start R.TPM.gene_end],1E6,chr_bin_lengths(1:22));

%% EXTRA


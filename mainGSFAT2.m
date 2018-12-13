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
if 1 == isTest
    load('testDataInfo.mat')
else
%     [dataInfo] = gsfatLoadUserInput;
    [dataInfo] = gsfatLoadUserInputV2;
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

%% save
save(sprintf('%s.mat',dataInfo.prjName),'dataInfo','H','R','-v7.3')

%% Analysis
%% Chromatin partitioning
[H] = gsfatChromPartition(dataInfo,H,R);

%% RNA-seq differential expression (between each sample)
% automatically selects the first and last time point, and between samples

%%%%%%%%% NEED TO FIX FROM HERE DOWN 12/8/18

% find samples to compare
samps = unique(dataInfo.sampleInfo.sample);



% run sample comparison
for iSampComp = 1:length(samps2comp)
    % create geneTable for DiffEx
    tempTreated = find(strncmp(R.TPM.Properties.VariableNames,'SW480_T72',9));
    tempUntreated = find(strncmp(R.TPM.Properties.VariableNames,'SW480_T0',8));
    
    meanTreated = mean(R.TPM{:,tempTreated},2);
    meanUntreated = mean(R.TPM{:,tempUntreated},2);
    meanBase = (meanTreated + meanUntreated) / 2;
    foldChange = meanTreated ./ meanUntreated;
    log2FC = log2(foldChange);
    
    R.geneTable = table(meanBase,meanTreated,meanUntreated,foldChange,log2FC);
    R.geneTable = [R.TPM(:,1:6),R.geneTable];
    
    [pvalue,padj] = matlabNegbinDE(R.expected_count{:,g1}, R.expected_count{:,g2});
    
    % add to the existing table
    R.geneTable.pvalue = pvalue;
    R.geneTable.padj = padj;
    
    % GSEA table %SRedit, 12/7/18
    % Sort genes based on padj, most increased to most decreased
    % create a temporary table to manipulate and organize
    tempTable = R.geneTable;
    tempTable(isnan(tempTable.log2FC),:) = [];
    tempTable(tempTable.padj >= 1,:) = [];
    tempTable(isnan(tempTable.padj),:) = [];
    
    % sort genes by padj, increased expression
    tempGseaTable = [];
    tempTable = sortrows(tempTable,'padj','ascend');
    tempGseaTable = tempTable(tempTable.log2FC > 0,:);
    
    % sort genes by padj, decreased expression
    tempTable = sortrows(tempTable,'padj','descend');
    tempGseaTable = [tempGseaTable;tempTable(tempTable.log2FC < 0,:)];
    
    % format table for GSEA
    % create "rnk" list. eg: [3 2 1 -1 -2 -3]
    rnkVal = sum(tempTable.log2FC > 0):-1:-sum(tempTable.log2FC < 0);
    rnkVal(rnkVal==0) = []; rnkVal = rnkVal';
    tempGseaTable = [tempGseaTable.geneName,table(rnkVal)];
    
    % write table to RNK file for GSEAQ input
    FileName = 'tcf7l2RnaseqGsea.rnk';
    writetable(tempGseaTable,'tcf7l2RnaseqGsea.rnk','Delimiter','\t',...
        'WriteVariableNames',1,'FileType','text')
    
    % format - add "#" at start of file
    S = fileread(FileName);
    S = ['#', S];
    FID = fopen(FileName, 'w');
    if FID == -1, error('Cannot open file %s', FileName); end
    fwrite(FID, S, 'char');
    fclose(FID);
    
end

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


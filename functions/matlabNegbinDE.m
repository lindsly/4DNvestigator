function [pvalue,padj] = matlabNegbinDE(sampleCounts1,sampleCounts2,varLink,lowCountThreshold,plotFlag)
%matlabNegbinDE calculates differential expression using a method similar
%to DESeq
%   This function is adapted from publicly available MATLAB scripts.
%   
%   Input
%   sampleCounts1:      column vector of raw counts for sample A
%   sampleCounts2:      column vector of raw counts for sample B
%   varLink:            link between variance for nbintest regression
%                       (Default: 'LocalRegression')
%   lowCountThreshold:  threshold for low number of counts (to be removed from model)
%   plotFlag:           indicator for output plots
%   
%   Output
%   pvalue:             p-values
%   padj:               adjusted p-values (Benjamini and Hochberg)
%
%   Scott Ronquist, scotronq@umich.edu. 3/7/19

%% Identifying Differentially Expressed Genes from RNA-Seq Data
% https://www.mathworks.com/help/bioinfo/examples/identifying-differentially...
% -expressed-genes-from-rna-seq-data.html
% This example shows how to test RNA-Seq data for differentially expressed
% genes using a negative binomial model.

% Copyright 2010-2016 The MathWorks, Inc.

%% default parameters
if nargin<3; varLink='LocalRegression'; end
if nargin<4; lowCountThreshold=10; end
if nargin<5; plotFlag=0; end

countAll = [sampleCounts1, sampleCounts2];

%% Inferring Differential Expression with a Negative Bionomial Model
% neg bin Regression
tRegression = nbintest(sampleCounts1,sampleCounts2,'VarianceLink',varLink);

%% p-value plots
if plotFlag
    h = plotVarianceLink(tRegression);
    
    % set custom title
    h(1).Title.String = 'Variance Link on Treated Samples';
    h(2).Title.String = 'Variance Link on Untreated Samples';
    
    % comparison
    h = plotVarianceLink(tRegression,'compare',true);
    
    % set custom title
    h(1).Title.String = 'Variance Link on Treated Samples';
    h(2).Title.String = 'Variance Link on Untreated Samples';
    
    % histogram of P-values
    figure;
    histogram(tRegression.pValue,100)
    xlabel('P-value')
    ylabel('Frequency')
    title('P-value enrichment')
end

%% Multiple Testing and Adjusted P-values
% compute the adjusted P-values (BH correction)
padj = mafdr(tRegression.pValue,'BHFDR',true);
pvalue = tRegression.pValue;

%% Filter low counts
lowCountGenes = all(countAll < lowCountThreshold, 2);

% plot p-value histogram
if plotFlag
    histogram(tRegression.pValue(~lowCountGenes),100)
    xlabel('P-value')
    ylabel('Frequency')
    title('P-value enrichment without low count genes')
end

padj(lowCountGenes) = NaN;
pvalue(lowCountGenes) = NaN;

end


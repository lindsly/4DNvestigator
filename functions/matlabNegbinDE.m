function [pvalue,padj] = matlabNegbinDE(sample_counts1,sample_counts2,var_link,lowCountThreshold)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%% Identifying Differentially Expressed Genes from RNA-Seq Data
% https://www.mathworks.com/help/bioinfo/examples/identifying-differentially...
% -expressed-genes-from-rna-seq-data.html
% This example shows how to test RNA-Seq data for differentially expressed
% genes using a negative binomial model.

% Copyright 2010-2016 The MathWorks, Inc.

%% default parameters
if nargin<4; lowCountThreshold=10; end
if nargin<3; var_link='LocalRegression'; end

count_all = [sample_counts1, sample_counts2];

%% Inferring Differential Expression with a Negative Bionomial Model
% neg bin Regression
t_regression = nbintest(sample_counts1,sample_counts2,'VarianceLink',var_link);
h = plotVarianceLink(t_regression);

% set custom title
h(1).Title.String = 'Variance Link on Treated Samples';
h(2).Title.String = 'Variance Link on Untreated Samples';

% comparison
h = plotVarianceLink(t_regression,'compare',true);

% set custom title
h(1).Title.String = 'Variance Link on Treated Samples';
h(2).Title.String = 'Variance Link on Untreated Samples';
  
%% histogram of P-values
figure;
histogram(t_regression.pValue,100)
xlabel('P-value')
ylabel('Frequency')
title('P-value enrichment')

% filter low counts
lowCountGenes = all(count_all{:,:} < lowCountThreshold, 2);
histogram(t_regression.pValue(~lowCountGenes),100)
xlabel('P-value')
ylabel('Frequency')
title('P-value enrichment without low count genes')

%% Multiple Testing and Adjusted P-values
% compute the adjusted P-values (BH correction)
padj = mafdr(t_regression.pValue,'BHFDR',true);
pvalue = t_regression.pValue;

end


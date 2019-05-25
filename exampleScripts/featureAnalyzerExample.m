%% 4DN Feature Analyzer example
% This example shows how the "4DN feature analyzer" can be used to find
% genes which change significantly in both structure and function over time
%
%   link to paper: (In preparation)
%
%   Version 1.0 (5/24/19)
%   Written by: Scott Ronquist
%   Contact:    scotronq@umich.edu
%   Created:    5/24/19
%
%   Revision History:
%   v1.0 (5/24/19)
%   * featureAnalyzerExample.m created

%% Load Data
% 4DNvestigator formatted Hi-C and RNA-seq data is available here:
%   https://drive.google.com/open?id=11XQ0CudxRvM5P6aI57yes2h8XXLhxlRr
clear
close all

load('myodData.mat')

%% 4DN Feature Analyzer Example
% select 4DN Feature Analyzer parameters
% This example is set up to analyze an entire genome at 1 Mb or 100 kb
% resolution
chrSelect = 11;
dimReduc = 'pca';
binSizeSelect = 1E5;

% extract region to analyze from Hi-C and RNA-seq
switch binSizeSelect
    case 1E6
        goiH = H.s1mb.oeTrim(H.s1mb.chrStartTrim(chrSelect):H.s1mb.chrStartTrim(chrSelect+1)-1,...
            H.s1mb.chrStartTrim(chrSelect):H.s1mb.chrStartTrim(chrSelect+1)-1,:);
        goiR = R.s1mb.tpmMeanTrim(H.s1mb.chrStartTrim(chrSelect):H.s1mb.chrStartTrim(chrSelect+1)-1,:);
        goi = R.s1mb.geneTrim(H.s1mb.chrStartTrim(chrSelect):H.s1mb.chrStartTrim(chrSelect+1)-1);
    case 1E5
        goiH = H.s100kb.oeTrim{chrSelect};
        goiR = R.s100kb.tpmMeanTrim{chrSelect};
        goi = R.s100kb.geneTrim{chrSelect};
    otherwise
        error('please select a correct bin size (1E6 or 1E5)')
end

% run the 4DNfeature analyzer
[features,score] = sfAnalysis(goiH,log2(goiR+1),goi,[],[],[],dimReduc);

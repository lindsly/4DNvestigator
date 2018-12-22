% This script performs all default analysis methods for the 4DNvestigator
% https://github.com/scotronq/4DNvestigator
%
% Variables
% dataInfo: MATLAB Structure that stores information that explains the input data.
% Includes data dependencies and file locations
%
% R: MATLAB Structure that contains all RNA-seq data
% H: MATLAB Structure that contains all Hi-C data
%
% sampleData is included wihtin the "sampleData" folder
%
% Scott Ronquist, 12/19/18

%% default start
clear
close all
restoredefaultpath
addpath(genpath('.'))
rmpath(genpath('./functions/depreciated'))

isTest = 0;

%% Input data Information
if isTest == 1
    load('testDataInfo.mat')
else
    [dataInfo] = fdnLoadUserInput;
end

%% load Hi-C
% loads and formats Hi-C data output from Juicer
if isTest == 1
    load('testDataInfo.mat','H')
else
    H = fdnLoadHic(dataInfo);
end

%% load RNA-seq
% loads and formats RNA-seq data output from RSEM
% H variable is loaded first and passed into this function to inform on
% which regions were "trimmed" (removed due to sparse coverage)
if isTest == 1
    load('testDataInfo.mat','R')
else
    R = fdnLoadRnaseq(dataInfo,H);
end

%% save
save(sprintf('./sampleData/data/%s.mat',dataInfo.projName),'dataInfo','H','R','-v7.3')

%% Chromatin partitioning
[H] = fdnChromPartition(dataInfo,H,R);

% find AB switch regions differences

%% RNA-seq differential expression (between each sample)
% automatically compares between all samples
[R] = fdnDiffExpGsaa(dataInfo,R);

%% 4DN feature analyzer - 1Mb level
% select sample
[features,score] = fdnSfAnalysis(dataInfo,R,H);

%% Matrix comparisons
% LP - genome-wide

% LP - chr by chr

% VNGE


%% EXTRA


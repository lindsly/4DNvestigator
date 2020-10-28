%% 4DN Feature Analyzer example
% This example script shows how each of the core functionalities of the
% 4DNvestigator can be called using default parameters
%
%   link to paper: (In preparation)
%
%   Written by: Scott Ronquist, Stephen Lindsly
%   Contact:    scotronq@umich.edu, lindsly@umich.edu

close all
clear
restoredefaultpath
addpath(genpath(pwd))

%% Load Data Example through the 4DNvestigator functions
if exist('myodData.mat','file')~=2 || exist('sampleMyodDataIndexTp-48_8_80.xlsx','file')~=2
    error(['Download all Hi-C and RNA-seq data ',...
        '<a href="https://drive.google.com/drive/folders/1xVjX7yqiOIPV_IfVKVDJJ79Ee0xMHGr8?usp=sharing">here</a>',...
        ' and save to the directory data\exampleData\'])
end

Index_Loc = 'data\exampleData\myodData\sampleMyodDataIndexTp-48_8_80.xlsx';
Data_Loc = 'data\exampleData\myodData\';

if exist('sampleMyodDataIndexTp-48_8_80.xlsx','file')==2
    [dataInfo] = fdnLoadUserInput(Index_Loc);
    [H] = fdnLoadHic(Data_Loc,dataInfo,'single');
    [R] = fdnLoadRnaseq(Data_Loc,dataInfo,H);
end

%% Feature Analyzer Example
% close all
Data_Loc = 'data\exampleData\myodData\myodData.mat';
Folder_Result = 'results\featureAnalyzerResults';
chrSelect = 11;
dimReduc = 'pca';
binSize = 1E5;

featureAnalyzerExample(Data_Loc, Folder_Result, chrSelect, dimReduc, binSize)

%% Simple Von Neumann Entropy Example
close all
Data_Loc = {'data\exampleData\4DNFIFLJLIS5.hic',...
            'data\exampleData\4DNFIOX3BGNE.hic'};
Folder_Result = 'results\vneExampleResults';
chrSelect = 14;
bpFrag = 'BP';
binSize = 1E5;

vneExample(Data_Loc, Folder_Result, chrSelect, bpFrag, binSize)

%% Expanded Von Neumann Entropy Example
close all
Data_Loc = {'data\exampleData\4DNFIFLJLIS5.hic';...
            'data\exampleData\4DNFIOX3BGNE.hic';...
            'https://hicfiles.s3.amazonaws.com/hiseq/imr90/in-situ/combined.hic';...
            'https://hicfiles.s3.amazonaws.com/hiseq/huvec/in-situ/combined.hic';...
            'https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic';...
            'http://hicfiles.s3.amazonaws.com/hiseq/rpe1/DarrowHuntley-2015/WT-combined.hic'};
Folder_Result = 'results\vneExampleExpandedResults';
chrSelect = 14;
bpFrag = 'BP';
binSize = 1E5;

vneExampleExpanded(Data_Loc, Folder_Result, chrSelect, bpFrag, binSize)

%% Larntz-Perlman Example
close all
Data_Loc = 'data\exampleData\myodData\myodData.mat';
Folder_Result = 'results\lpExampleResults';

lpExample(Data_Loc, Folder_Result);

%% Network Example
close all
Data_Loc = 'data\exampleData\networkData\'; %%% data input folder
Folder_Result = 'results\networkExamplesResults'; %%% output result file name

networkExamples(Data_Loc, Folder_Result);

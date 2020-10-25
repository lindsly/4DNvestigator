close all, clear all
restoredefaultpath
addpath(genpath(pwd))

% TEMPORARY
addpath(genpath('\\172.17.109.24\internal_4dn\projects\4DNvestigator_data'))

%% Load Data Example through the 4DNvestigator functions
if exist('myodData.mat','file')~=2 || exist('sampleMyodDataIndexTp-48_8_80.xlsx','file')~=2
    error(['Download time series Hi-C and RNA-seq data ',...
        '<a href="https://drive.google.com/open?id=1lSyU-7I0ME3X70Mt_-HjLMtPc-BMKHxm">here</a>'])
end

Index_Loc = '\\172.17.109.24\internal_4dn\projects\4DNvestigator_data\myodData\sampleMyodDataIndexTp-48_8_80.xlsx';
Data_Loc = '\\172.17.109.24\internal_4dn\projects\4DNvestigator_data\myodData\';


if exist('sampleMyodDataIndexTp-48_8_80.xlsx','file')==2
    [dataInfo] = fdnLoadUserInput(Index_Loc);
    [H] = fdnLoadHic(Data_Loc,dataInfo,'single');
    [R] = fdnLoadRnaseq(Data_Loc,dataInfo,H);
end

%% Feature Analyzer Example
close all
% *v* Temporary storage solution *v*
Data_Loc = '\\172.17.109.24\internal_4dn\projects\4DNvestigator_data\myodData\myodData.mat';
Folder_Result = 'results\featureAnalyzerResults';
chrSelect = 11;
dimReduc = 'pca';
binSize = 1E5;
featureAnalyzerExample(Data_Loc, Folder_Result, chrSelect, dimReduc, binSize)

%% Simple Von Neumann Entropy Example
close all
Data_Loc = {'\\172.17.109.24\internal_4dn\projects\4DNvestigator_data\4DNFIFLJLIS5.hic',...
            '\\172.17.109.24\internal_4dn\projects\4DNvestigator_data\4DNFIOX3BGNE.hic'};
Folder_Result = 'results\vneExampleResults';
chrSelect = 14;
bpFrag = 'BP';
binSize = 1E5;

vneExample(Data_Loc, Folder_Result, chrSelect, bpFrag, binSize)

%% Expanded Von Neumann Entropy Example
close all
Data_Loc = {'\\172.17.109.24\internal_4dn\projects\4DNvestigator_data\4DNFIFLJLIS5.hic';...
            '\\172.17.109.24\internal_4dn\projects\4DNvestigator_data\4DNFIOX3BGNE.hic';...
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
Data_Loc = '\\172.17.109.24\internal_4dn\projects\4DNvestigator_data\myodData\myodData.mat';
Folder_Result = 'results\lpExampleResults';

lpExample(Data_Loc, Folder_Result);

%% Network Example
close all
% The folder "networkData" and its contents must be downloaded to the current
% directory to run this function. "networkData" can be downloaded here:
% https://drive.google.com/drive/folders/17XhveY8HDjeh3KWz43BBlCu1xN8yVs16?usp=sharing

Data_Loc = '\\172.17.109.24\internal_4dn\projects\4DNvestigator_data\networkData\'; %%% data input folder
Folder_Result = 'results\networkExamplesResults'; %%% output result file name

networkExamples(Data_Loc, Folder_Result);
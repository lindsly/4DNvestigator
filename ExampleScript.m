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

if exist('sampleMyodDataIndexTp-48_8_80.xlsx','file')==2
    [dataInfo] = fdnLoadUserInput(indexFile);
    [H] = fdnLoadHic(dataInfo,'single');
    [R] = fdnLoadRnaseq(dataInfo,H);
end

%% Feature Analyzer Example
% *v* Temporary storage solution *v*
Folder_Data = '\\172.17.109.24\internal_4dn\projects\4DNvestigator_data\myodData';
% Folder_Data = 'exampleData\myodData';
Folder_Result = 'results\featureAnalyzerResults';
chrSelect = 11;
dimReduc = 'pca';
binSize = 1E5;
featureAnalyzerExample(Folder_Data, Folder_Result, chrSelect, dimReduc, binSize)

%% Simple Von Neumann Entropy Example

%% Extended Von Neumann Entropy Example

%% Larntz-Perlman Example

%% Network Example
% The folder "networkData" and its contents must be downloaded to the current
% directory to run this function. "networkData" can be downloaded here:
% https://drive.google.com/drive/folders/17XhveY8HDjeh3KWz43BBlCu1xN8yVs16?usp=sharing

Folder_Data = sprintf('networkData');  %%%   CellTrans_Data data file name
Folder_Result = sprintf('Figs_J_bioinfo'); %%% output result file name

networkExamples(Folder_Data, Folder_Result);
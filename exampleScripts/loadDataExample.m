%% Description
% This script provides and overview of how the 4DNvestigator loads data
%
% Scott Ronquist, scotronq@umich.edu. 4/29/19

%% Set up
clear
close all

% Add 4DNvestigator tools to path
filepath = mfilename('fullpath');
fdnPath = filepath(1:strfind(filepath,'4DNvestigator')+12);
addpath(genpath(fdnPath))

%% Select Data set to Load
% Time Series Hi-C and RNA-seq data available at the following link:
%   https://drive.google.com/open?id=1lSyU-7I0ME3X70Mt_-HjLMtPc-BMKHxm

if exist('myodData.mat','file')==2 && exist('sampleMyodDataIndexTp-48_8_80.xlsx','file')==2
    error(['Download time series Hi-C and RNA-seq data ',...
        '<a href="https://drive.google.com/open?id=1lSyU-7I0ME3X70Mt_-HjLMtPc-BMKHxm">here</a>'])
end

%% Load data through the 4DNvestigator functions
if exist('sampleMyodDataIndexTp-48_8_80.xlsx','file')==2
    [dataInfo] = fdnLoadUserInput(indexFile);
    [H] = fdnLoadHic(dataInfo,'single');
    [R] = fdnLoadRnaseq(dataInfo,H);
end
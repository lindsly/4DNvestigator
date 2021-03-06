%% 4DN example
% This example script shows how each of the core functionalities of the
% 4DNvestigator can be called using default parameters
%
%   link to paper: (In preparation)
%
%   Written by: Stephen Lindsly, Scott Ronquist 
%   Contact:    lindsly@umich.edu, scotronq@umich.edu

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

%% Load Generic Hi-C data (pre-processed chromosomes)
numChr = 24;
H_test = fdnEmptyHMatrix(numChr);

Data_Loc = 'data\exampleData\GenericHiCData\';
File_Name = 'generic_chr1_hic_test.txt';
binSize = 1E5;
chrSelect = 1;
H_test = fdnLoadGenericHic(Data_Loc, File_Name, chrSelect, binSize, 1, H_test);

%% Feature Analyzer Example
% close all
Data_Loc = 'data\exampleData\myodData\myodData.mat';
Folder_Result = 'results\featureAnalyzerResults';
chrSelect = 11;
dimReduc = 'tsne';
binSize = 1E5;
sfType = 'sfmatrix'; % 'phaseplane' can replace 'sfmatrix' to plot a phase plane
featureType = 'rna'; % 'other' can replace 'rna' for other 1D genomic features

featureAnalyzerExample(Data_Loc, Folder_Result, chrSelect, dimReduc, binSize, sfType, featureType)

%% Simple Von Neumann Entropy Example
close all
Data_Loc = {'data\exampleData\4DNFIFLJLIS5.hic',...
            'data\exampleData\4DNFIOX3BGNE.hic'};
Folder_Result = 'results\entropyExampleResults';
chrSelect = 14;
bpFrag = 'BP';
binSize = 1E5;

entropyExample(Data_Loc, Folder_Result, chrSelect, bpFrag, binSize)

%% Expanded Von Neumann Entropy Example
close all
Data_Loc = {'data\exampleData\4DNFIFLJLIS5.hic';...
            'data\exampleData\4DNFIOX3BGNE.hic';...
            'https://hicfiles.s3.amazonaws.com/hiseq/imr90/in-situ/combined.hic';...
            'https://hicfiles.s3.amazonaws.com/hiseq/huvec/in-situ/combined.hic';...
            'https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic';...
            'http://hicfiles.s3.amazonaws.com/hiseq/rpe1/DarrowHuntley-2015/WT-combined.hic'};
Folder_Result = 'results\entropyExampleExpandedResults';
chrSelect = 14;
bpFrag = 'BP';
binSize = 1E5;

entropyExampleExpanded(Data_Loc, Folder_Result, chrSelect, bpFrag, binSize)

%% Larntz-Perlman Example
close all
Data_Loc = 'data\exampleData\myodData\myodData.mat';
Folder_Result = 'results\lpExampleResults';

lpExample(Data_Loc, Folder_Result);

%% Tensor Entropy Example
% close all
Data_Loc = 'data\exampleData\tensorData\myodEntropy.mat';
Folder_Result = 'results\tensorEntropyExampleResults';

tensorEntropyExample(Data_Loc, Folder_Result);

function [dataInfo] = fdnLoadUserInput(indexFile,projectName,outputPath)
%fdnLoadUserInput loads User input data
%   This function loads formats data structures for 4DNvestigator analysis
%
%   Input
%   indexFile:      file that contains information on all samples to be
%                   processed (string)
%   projectName:    name of project
%   outputPath:     folder for project analysis output
%
%   Output
%   dataInfo:       structure that contains all relevant experiment data information
%
%   Version 1.1 (06/11/19)
%   Written by: Scott Ronquist
%   Contact:    scotronq@umich.edu
%   Created:    1/22/18
%   
%   Revision History:
%   v1.0 (1/22/18)
%   v1.1 (1/22/18)
%   * removed projectName and outputPath sections (now in fdnSave.m)

%% Get computer info
dataInfo.delim = filesep;
dataInfo.hicCMap = [ones(64,1),[1:-1/63:0]',[1:-1/63:0]'];

%% Get Project Name
% if exist('projectName','var')
%     dataInfo.projName = projectName;
% else
%     temp = inputdlg('Input Project Name:','Project Name');
%     dataInfo.projName = temp{1};
% end

%% Create analysis output directory
% if exist('outputPath','var')
%     dataInfo.path.output = outputPath;
% else
%     fprintf('Select Output folder\n')
%     dataInfo.path.output = uigetdir(pwd,'Select Output folder');
% end
% 
% % make subdirectories
% mkdir([dataInfo.path.output,dataInfo.delim,'figures'])
% mkdir([dataInfo.path.output,dataInfo.delim,'tables'])
% mkdir([dataInfo.path.output,dataInfo.delim,'data'])
% mkdir([dataInfo.path.output,dataInfo.delim,'data',dataInfo.delim,'gsaa'])

%% Load from input index file
if exist('indexFile','var')
    if ischar(indexFile)
        dataInfo.indexFile = indexFile;
        
        % webread file if AWS
        if contains(indexFile,'http')
            dataInfo.sampleInfo = webread(dataInfo.indexFile);
        else
            dataInfo.sampleInfo = readtable(dataInfo.indexFile);
        end
        
    elseif istable(indexFile)
        dataInfo.sampleInfo = indexFile;
    end
else
    answer = [];
    while ~strcmp(answer,'Yes')
        answer = questdlg('Do you have an Index File for your data set?', ...
            'input Index File',...
            'Yes', 'No', 'What is this?',...
            'Yes');
        
        % Handle response
        switch answer
            case 'Yes'
                [file,path] = uigetfile('*','Select Index File:',pwd);
                dataInfo.indexFile = [path,file];
                
            case 'No'
                fprintf('please create one, manual input not available at this time.\n')
                
            case 'What is this?'
                fprintf(['An index file is a tab or comma separated file that includes ',...
                    'information on the samples being process.\n',...
                    'files include the following fields: path, dataType, sample, timePoint.\n\n',...
                    'path: full path for the data file on your computer.\n',...
                    'dataType: type of data (hic or rnaseq).\n',...
                    'sample: sample name (string).\n',...
                    'timePoint: time point of data (integer).\n\n',...
                    'An example Index File is located in ./sampleData/sampleDataIndex.xlsx\n'])
        end
    end
    % read index file
    dataInfo.sampleInfo = readtable(dataInfo.indexFile);
end

% extract sample names from sample filename
tempName = cell(height(dataInfo.sampleInfo),1);
for iSample = 1:height(dataInfo.sampleInfo)
    [filepath,name,ext] = fileparts(dataInfo.sampleInfo.path{iSample});
    nameSimple = strsplit(name,{'.'});
    tempName{iSample} = nameSimple{1};
end
dataInfo.sampleInfo.name = tempName;

%% Read .hic header files
fprintf('reading .hic header information...\n')

% get .hic samples locations
hicSampleLocs = find(ismember(dataInfo.sampleInfo.dataType,'hic'));

% read all .hic headers and get reference genome information
hicHeader = cell(length(hicSampleLocs),1);
refGenome = cell(length(hicSampleLocs),1);
idxRefGenome = [];
for iHicSample = 1:length(hicSampleLocs)
    if ismember('refGenome',dataInfo.sampleInfo.Properties.VariableNames)
        idxRefGenome = dataInfo.sampleInfo.refGenome{hicSampleLocs(iHicSample)};
    else
        idxRefGenome = [];
    end
    
    % read .hic header
    hicHeader{iHicSample} = readHicHeader(dataInfo.sampleInfo.path{hicSampleLocs(iHicSample)},...
        idxRefGenome);
    refGenome{iHicSample} = hicHeader{iHicSample}.refGenome;
end

% determine if all samples were mapped to same reference genome
if length(unique(refGenome)) == 1
    dataInfo.hicHeader = hicHeader{iHicSample};
    
    % remove "ALL" and "M" chr
    dataInfo.hicHeader.Chromosomes(ismember(upper(dataInfo.hicHeader.Chromosomes.chr),...
        {'ALL','M','MT'}),:) = [];
    
else
    error('.hic files have different reference genomes, cannot compare')
end

%% Create replicate identifier
[uvals, ~, uidx] = unique(strcat(dataInfo.sampleInfo.sample,...
    cellstr(num2str(dataInfo.sampleInfo.timePoint)), dataInfo.sampleInfo.dataType));
replicate = ones(height(dataInfo.sampleInfo),1);    % used to mostly to copy the class and size
for K = 1 : length(uvals)
  mask = uidx == K;
  replicate(mask) = [1:sum(mask)]';
end
dataInfo.sampleInfo.replicate = replicate;

% Create uniqname from info
dataInfo.sampleInfo.uniqueName = strcat(dataInfo.sampleInfo.dataType,...
    {'_s'},dataInfo.sampleInfo.sample,...
    {'_t'},cellstr(num2str(dataInfo.sampleInfo.timePoint)),...
    {'_r'},cellstr(num2str(dataInfo.sampleInfo.replicate)));

% sort sampleInfo
dataInfo.sampleInfo = sortrows(dataInfo.sampleInfo,...
    {'dataType','sample','timePoint'},{'ascend','ascend','ascend'});

% Create sampleSlice (to be filled later as data is loaded)
dataInfo.sampleInfo.index = [1:sum(ismember(dataInfo.sampleInfo.dataType,'hic')),...
    1:sum(ismember(dataInfo.sampleInfo.dataType,'rnaseq'))]';

end


function [dataInfo] = gsfatLoadUserInputV2()
%gsfatUserInput loads User input data
%   This function load user input data for downstream analysis

%% get computer info
if isunix
    dataInfo.delim = '/';
else
    dataInfo.delim = '\';
end


%% load from input index file
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
            dataInfo.sampleInfo = readtable(dataInfo.indexFile);
            
            % extract sample name from sample filename
            tempName = cell(height(dataInfo.sampleInfo),1);
            for iSample = 1:height(dataInfo.sampleInfo)
                temp1 = strsplit(dataInfo.sampleInfo.path{iSample},dataInfo.delim);
                temp2 = strsplit(temp1{end},{'.'});
                tempName{iSample} = temp2{1};
            end
            dataInfo.sampleInfo.name = tempName;
            
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

%% get available information from .hic files
% refgenome
fprintf('reading .hic header information...\n')
hicSampleLocs = find(ismember(dataInfo.sampleInfo.dataType,'hic'));

hicHeader = cell(length(hicSampleLocs),1);
refGenome = cell(length(hicSampleLocs),1);
for iHicSample = 1:length(hicSampleLocs)
    hicHeader{iHicSample} = readHicHeader(dataInfo.sampleInfo.path{hicSampleLocs(iHicSample)});
    refGenome{iHicSample} = hicHeader{iHicSample}.refGenome;
end

if length(unique(refGenome)) == 1
    dataInfo.hicHeader = hicHeader{iHicSample};
    dataInfo.refGenome = dataInfo.hicHeader.refGenome;
    
    % need to fix this later vvvvvvvvvv
    dataInfo.numChr = 22;
    dataInfo.chrSizes = readtable(sprintf('%s.chrom.sizes',dataInfo.refGenome),'filetype','text');
    dataInfo.chrSizes = dataInfo.chrSizes(1:dataInfo.numChr,:);
    % need to fix this later ^^^^^^^^^^
else
    error('.hic files have different reference genomes')
end

%% get peripheral info
% Get experiment Name
temp = inputdlg('Input Project Name:','Project Name');
dataInfo.projName = temp{1};

% Create replicate identifier
[uvals, ~, uidx] = unique(strcat(dataInfo.sampleInfo.sample,dataInfo.sampleInfo.dataType));
replicate = ones(height(dataInfo.sampleInfo),1);    % used to mostly to copy the class and size
for K = 1 : length(uvals)
  mask = uidx == K;
  replicate(mask) = [1:sum(mask)]';
end
dataInfo.sampleInfo.replicate = replicate;

% Create uniqname from info
dataInfo.sampleInfo.uniqueName = strcat(dataInfo.sampleInfo.dataType,...
    {'_s'},dataInfo.sampleInfo.sample,{'_r'},cellstr(num2str(replicate)));

%% Format output directory
% get output directory
dataInfo.path.output = uigetdir(pwd,'Select Output folder');
mkdir([dataInfo.path.output,dataInfo.delim,'figures'])
mkdir([dataInfo.path.output,dataInfo.delim,'/tables'])

%% depreciated
% % select reference genome
% temp = {'hg19','hg38'};
% [indx,~] = listdlg('PromptString','Select a Reference Genome:',...
%     'SelectionMode','single','ListString',temp);
% dataInfo.refGenome = temp{indx};
% 
% 
% dataInfo.numChr = 22;
% dataInfo.chrSizes = readtable(sprintf('%s.chrom.sizes',dataInfo.refGenome),'filetype','text');
% dataInfo.chrSizes = dataInfo.chrSizes(1:dataInfo.numChr,:);

end


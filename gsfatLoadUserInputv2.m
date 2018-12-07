function [dataInfo] = gsfatLoadUserInputv2()
%gsfatUserInput loads User input data
%   This function load user input data for downstream analysis

%% load from input file
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
            dataInfo.sampleInfo.name = file;
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

%% get periforal info
% select reference genome
temp = {'hg19','hg38'};
[indx,~] = listdlg('PromptString','Select a Reference Genome:',...
    'SelectionMode','single','ListString',temp);
dataInfo.refGenome = temp{indx};

% Experiment Name
temp = inputdlg('Input Project Name:','Project Name');
dataInfo.projName = temp{1};

% create replicate identifier
[uvals, ~, uidx] = unique(strcat(unique(dataInfo.sampleInfo.sample),dataInfo.sampleInfo.dataType));
replicate = ones(height(dataInfo.sampleInfo),1);    % used to mostly to copy the class and size
for K = 1 : length(uvals)
  mask = uidx == K;
  replicate(mask) = [1:sum(mask)]';
end

% create uniqname from info
uniqname = strcat(dataInfo.sampleInfo.dataType,{'_s'},sample,{'_r'},cellstr(num2str(replicate)));
dataInfo.sampleInfo = [dataInfo.sampleInfo,table(sample),table(replicate),table(uniqname)];



numChr = 22;

%% Input data Information
% select your input data folders
temp = inputdlg('Input Project Name:','Project Name');
dataInfo.projName = temp{1};
dataInfo.path.rnaseq = uigetdir(pwd,'Select RNA-seq data folder:');
dataInfo.path.hic = uigetdir(pwd,'Select Hi-C data folder:');
dataInfo.path.output = uigetdir(pwd,'Select output directory:');

% select reference genome
temp = {'hg19','hg38'};
[indx,tf] = listdlg('PromptString','Select a Reference Genome:',...
    'SelectionMode','single','ListString',temp);
dataInfo.refGenome = temp{indx};

% select Hi-C Processing programs
temp = {'juicer'};
[indx,tf] = listdlg('PromptString','Select a Hi-C processing program:',...
    'SelectionMode','single','ListString',temp);
dataInfo.dataFormat.hic = temp{indx};

% select Hi-C and RNa-seq Processing programs - can be replaced by filetype
% extension in the future
temp = {'RSEM'};
[indx,tf] = listdlg('PromptString','Select a Hi-C processing program:',...
    'SelectionMode','single','ListString',temp);
dataInfo.dataFormat.rnaseq = temp{indx};

% add to path data folders
addpath(genpath(dataInfo.path.rnaseq))
addpath(genpath(dataInfo.path.hic))

%% organize samples
% create sample info table, define samples/replicates
listing=dir(dataInfo.path.hic);[listing(:).dataType] = deal('hic');
dataInfo.sampleInfo=struct2table(listing(3:end));

listing=dir(dataInfo.path.rnaseq);[listing(:).dataType] = deal('rnaseq');
dataInfo.sampleInfo=[dataInfo.sampleInfo;struct2table(listing(3:end))];
sample = inputdlg(dataInfo.sampleInfo.name,'define samples',[1 50]);

% create replicate identifier
[uvals, ~, uidx] = unique(strcat(sample,dataInfo.sampleInfo.dataType));
replicate = ones(height(dataInfo.sampleInfo),1);    % used to mostly to copy the class and size
for K = 1 : length(uvals)
  mask = uidx == K;
  replicate(mask) = [1:sum(mask)]';
end

% create uniqname from info
uniqname = strcat(dataInfo.sampleInfo.dataType,{'_s'},sample,{'_r'},cellstr(num2str(replicate)));
dataInfo.sampleInfo = [dataInfo.sampleInfo,table(sample),table(replicate),table(uniqname)];

% get chromosome sizes
dataInfo.chrSizes = readtable(sprintf('%s.chrom.sizes',dataInfo.refGenome),'filetype','text');
dataInfo.chrSizes = dataInfo.chrSizes(1:numChr,:);

%% format output directory
mkdir([dataInfo.path.output,'/figures'])
mkdir([dataInfo.path.output,'/tables'])

end


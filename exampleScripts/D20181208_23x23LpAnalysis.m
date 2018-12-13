% Goal: perform Larntz Perlman analysis for 23 x 23 correlation matrix:

%% load data/paths
clear
close all
addpath(genpath('E:\MATLAB\gsfat'))

%% select cell types and Hi-C parameters
%%% PARAMETERS vvv
param.binType = 'BP';
param.binSize = 1E6;
param.norm1d = 'KR';
param.norm3d = 'oe';
param.intraFlag = 1;
param.chr = 14;
param.sampleSelect = {'IMR90','HFFc6','H1-hESC'};
%%% PARAMETERS ^^^

fn = {'http://hicfiles.s3.amazonaws.com/hiseq/rpe1/DarrowHuntley-2015/WT-combined.hic';...
    'https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic';...
    'https://hicfiles.s3.amazonaws.com/hiseq/imr90/in-situ/combined.hic';...
    'https://hicfiles.s3.amazonaws.com/hiseq/hmec/in-situ/combined.hic';...
    'https://hicfiles.s3.amazonaws.com/hiseq/nhek/in-situ/combined.hic';...
    'https://hicfiles.s3.amazonaws.com/hiseq/huvec/in-situ/combined.hic';...
    'https://data.4dnucleome.org/files-processed/4DNFIFLJLIS5/@@download/4DNFIFLJLIS5.hic';...
    'https://data.4dnucleome.org/files-processed/4DNFIOX3BGNE/@@download/4DNFIOX3BGNE.hic'};

cellType = {'RPE1 WT';'GM12878';'IMR90';'HMEC';'NHEK';'HUVEC';'HFFc6';'H1-hESC'};
dataInfo = table(cellType,fn);
dataInfo = sortrows(dataInfo,'cellType','ascend');
dataInfo = dataInfo(ismember(dataInfo.cellType,param.sampleSelect),:);

%% load hic
H = zeros(22,22,height(dataInfo));
for iSample = 1:height(dataInfo)
    for iChr1 = 1:22
        for iChr2 = 1:22
            disp([iSample,iChr1,iChr2])
            tempH = hic2mat('observed','none',dataInfo.fn{iSample},...
                iChr1,iChr2,param.binType,param.binSize,[]);
            H(iChr1,iChr1,iSample) = sum(tempH(:));
        end
    end
end

%% normalize data
% find chr start and end
chrStart = find(diff(sum(logical(H{iSample}))) > 50);

%% visualize data


%% LP analyze data





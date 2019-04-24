% Goal: load Hi-C data from our server onto alternate laptop 
clear
close all

%% load paths
restoredefaultpath
addpath(genpath('.'))

%% Load Hi-C data
hicLoc = '\\172.17.109.24\internal_4DN\projects\tcf7l2_silence_sw480_and_rpe\hic\processed\hg19\Sample_76099\aligned\inter_30.hic';
%hicLoc = 'https://drive.google.com/open?id=1p94KgbbQ1lA_xiFClA80ARINM5pH8h6A';

%%% PARAMETERS vvv
hicParam.binType = 'BP';
hicParam.binSize = 1E5;
hicParam.norm1d = 'KR';
hicParam.norm3d = 'oe';
hicParam.intraFlag = 1;
hicParam.chr = 1;
%%% PARAMETERS ^^^

[H,hicHeader] = hic2mat(hicParam.norm3d,hicParam.norm1d,hicLoc,...
    hicParam.chr,hicParam.chr,hicParam.binType,hicParam.binSize,hicParam.intraFlag,1);

figure, imagesc(log(H))

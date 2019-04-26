% Goal: This script performs the 4DN feature analyzer analysis on MYOD
% reprogramming time points
%
% Info
%   Sample_64585    control
%   Sample_71530    T1
%   Sample_71531    T12
%
% Scott Ronquist, scotronq@umich.edu. 4/23/19

%% set up
clear
close all

%% load data
if ~isfile('./data/myodTsData.mat')
    [dataInfo] = fdnLoadUserInput('myodDataIndex.xlsx','myod','.');
    [H] = fdnLoadHic(dataInfo);
    [R] = fdnLoadRnaseq(dataInfo,H);
    
    save('./data/myodTsData','H','R','dataInfo','-v7.3')
else
    load('./data/myodTsData.mat')
end

%% select ROI
% select 100kb, chr
chrSelect = 11;
goiH = H.s100kb.oeTrim{chrSelect};
goiR = R.s100kb.tpmMeanTrim{chrSelect};
goi = R.s100kb.geneTrim{chrSelect};

%% 4dn feature analyzer
% [features,score] = fdnSfAnalysis(dataInfo,R,H,sampleSelect,goi,hExtract);
dimReduc = 'pca';
[features,score] = sfAnalysis(goiH,goiR,goi,[],[],[],dimReduc);
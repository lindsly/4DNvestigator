% Goal: example of time series differential expression
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

%% time series differential expression

[R] = fdnDiffExpGsaa(dataInfo,R);
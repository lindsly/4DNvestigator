% Goal: test Larntz-Perlman method for comparing covariance matrices
%
%   Reference: Koziol, James A., et al. "A graphical technique for
%   displaying correlation matrices." The American Statistician 51.4
%   (1997): 301-304.

%% Load 
clear
close all

restoredefaultpath
addpath(genpath('.'))

%% simple test
% parameters
rng(1)
alphaParam = .95;
P = 5;
N = 15;
noise = 1;

dataA = rand(P,N);
dataB = dataA+rand(P,N)*noise;

corrA = corr(dataA');
corrB = corr(dataB');

R = cat(3,corrA,corrB);

[H0,M] = larntzPerlman(R,N,alphaParam);

% figure
figure('position',[100 100 900 900])
tempMin = min([min(corrA(:)),min(corrB(:))]);
subplot(2,2,1), imagesc(corrA), axis square, caxis([tempMin 1]), colorbar, title('Correlation of A')
subplot(2,2,2), imagesc(corrB), axis square, caxis([tempMin 1]), colorbar, title('Correlation of B')
subplot(2,2,3), imagesc(abs(corrA-corrB)), axis square,  colorbar, title('A - B')
subplot(2,2,4), imagesc(M), axis square, colorbar, title('Larntz-Perlman')%title(sprintf('%.4f'))

%% Hi-C type example
%%% PARAMETERS vvv
hicParam.binType = 'BP';
hicParam.binSize = 1E6;
hicParam.norm1d = 'KR';
hicParam.norm3d = 'oe';
hicParam.intraFlag = 1;
hicParam.chr = 14;
sampleSelect = {'IMR90','HFFc6'};
%%% PARAMETERS ^^^

sampleDataLoc = {'http://hicfiles.s3.amazonaws.com/hiseq/rpe1/DarrowHuntley-2015/WT-combined.hic',...
    'https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic',...
    'https://hicfiles.s3.amazonaws.com/hiseq/imr90/in-situ/combined.hic',...
    'https://hicfiles.s3.amazonaws.com/hiseq/hmec/in-situ/combined.hic',...
    'https://hicfiles.s3.amazonaws.com/hiseq/nhek/in-situ/combined.hic',...
    'https://hicfiles.s3.amazonaws.com/hiseq/huvec/in-situ/combined.hic',...
    'https://data.4dnucleome.org/files-processed/4DNFIFLJLIS5/@@download/4DNFIFLJLIS5.hic',...
    'https://data.4dnucleome.org/files-processed/4DNFIOX3BGNE/@@download/4DNFIOX3BGNE.hic'};

sampleNames = {'RPE1 WT','GM12878','IMR90','HMEC','NHEK','HUVEC','HFFc6','H1-hESC'};

sampleDataLoc = sampleDataLoc(ismember(sampleNames,sampleSelect));
sampleNames = sampleNames(ismember(sampleNames,sampleSelect));

H1 = hic2mat(hicParam.norm3d,hicParam.norm1d,sampleDataLoc{1},...
    hicParam.chr,hicParam.chr,hicParam.binType,hicParam.binSize,hicParam.intraFlag);
H2 = hic2mat(hicParam.norm3d,hicParam.norm1d,sampleDataLoc{2},...
    hicParam.chr,hicParam.chr,hicParam.binType,hicParam.binSize,hicParam.intraFlag);

[HtrimAll,badLocs] = hicTrim(padconcatenation_sr(H1,H2,3),1,.5);

H1corr = corr(HtrimAll(:,:,1));
H2corr = corr(HtrimAll(:,:,2));

% Larntz-Perlman
R = cat(3,H1corr,H2corr);
[H0,M,P] = larntzPerlman(R,size(R,1),alphaParam);

mask = triu(true(size(M)),1);
r = sort(M(mask),'descend');
% M = (M >= r(5));

% P-val thresh matrix
figure, imagesc(P < .01)


figure('position',[100 100 900 900])
tempMin = min([min(H1corr(:)),min(H2corr(:))]);
subplot(2,2,1), imagesc(H1corr), axis square, caxis([tempMin 1]), colorbar, title(sampleNames{1})
subplot(2,2,2), imagesc(H2corr), axis square, caxis([tempMin 1]), colorbar, title(sampleNames{2})
subplot(2,2,3), imagesc(abs(H1corr-H2corr)), axis square,  colorbar, title(sprintf('%s - %s',sampleNames{1},sampleNames{2}))
subplot(2,2,4), imagesc(M), axis square, colorbar, title('Larntz-Perlman')%title(sprintf('%.4f'))
linkaxes(get(gcf,'children'))

%% overview of method
% "A matrix of pair-wise correlations of the P300 amplitudes (i.e., the
% amount of strength of the electrophysiological response) at the 19 sites
% across the subjects was then calculated for each laboratory"

%% EXTRA
% %% create sample data
% % parameters
% rng(1)
% numSubjectsPerLab = 15;
% numLabs = 6;
% numElectrodes = 19;
% labVar = 1;
% subjectVar = 1;
% electrodeVar = 1;
% electrodeMean = rand(numElectrodes,1)*100;
% 
% % create sample data
% corrData = zeros(numElectrodes,numElectrodes,numLabs);
% 
% for iLab = 1:numLabs
%     labNoise = randn*labVar;
%     electrodeData = zeros(numElectrodes,numSubjectsPerLab);
%     
%     for iSubject = 1:numSubjectsPerLab
%         subjectNoise = randn*subjectVar;
%         
%         for iElectrode = 1:numElectrodes
%             electrodeData(iElectrode,iSubject) =...
%                 normrnd(electrodeMean(iElectrode),electrodeVar) + subjectNoise + labNoise;
%         end
%     end
%     corrData(:,:,iLab) = corr(electrodeData');%figure, imagesc(corr(electrodeData'))
% end
% %figure, imagesc(corrData(:,:,1))
% 
% %% test larntz-Perlman on sample data
% alphaParam = .95;
% [H0,z] = larntzPerlman(corrData, numSubjectsPerLab, alphaParam);
% 
% %% Hi-C data
% 
% 
% %% trisomy 7 data



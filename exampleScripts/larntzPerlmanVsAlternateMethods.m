%   This script compares how the Larntz-Perlman procedure performs relative
%   to alternate Hi-C comparison methods
%
%   References:
%   - Characterizing the 3D structure and dynamics of chromosomes and
%   proteins in a common contact matrix framework.
%   https://academic.oup.com/nar/article/46/16/8143/5051111
%
%   -HiCRep: assessing the reproducibility of Hi-C data using a stratum-
%   adjusted correlation coefficient.
%   https://genome.cshlp.org/content/early/2017/08/30/gr.220640.117
%
%   - GenomeDISCO: a concordance score for chromosome conformation capture
%   experiments using random walks on contact map graphs.
%   https://academic.oup.com/bioinformatics/article/34/16/2701/4938489
%
%   - HiC-spector: a matrix library for spectral and reproducibility
%   analysis of Hi-C contact maps.
%   https://academic.oup.com/bioinformatics/article/33/14/2199/3078603
%
%
%   Version 1.0 (4/14/19)
%   Written by: Scott Ronquist
%   Contact: scotronq@umich.edu
%   Contributors:
%   Created: 4/14/19
%   Revision History:

%% Script set-up
clear
close all

%% load data
% Public Hi-C data locations
names = {'GM12878 Primary','GM12878 Replicate','IMR90','NHEK','HUVEC'};
dataLoc = {'https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/primary.hic',...
    'https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/replicate.hic',...
    'https://hicfiles.s3.amazonaws.com/hiseq/imr90/in-situ/combined.hic',...
    'https://hicfiles.s3.amazonaws.com/hiseq/nhek/in-situ/combined.hic',...
    'https://hicfiles.s3.amazonaws.com/hiseq/huvec/in-situ/combined.hic'};
sampleInfo = table(names,dataLoc);

%% Hi-C type example
% Define Hi-C extraction parameters
%%% PARAMETERS vvv
hicParam.binType = 'BP';
hicParam.binSize = 1E6;
hicParam.norm1d = 'KR';
hicParam.norm3d = 'oe';
hicParam.intraFlag = 1;
hicParam.chr = 1;
%%% PARAMETERS ^^^

%% Load Hi-C
H = [];
for iSample = 1:length(sampleNames)
    tempH = hic2mat(hicParam.norm3d,hicParam.norm1d,dataLoc{iSample},...
        hicParam.chr,hicParam.chr,hicParam.binType,hicParam.binSize,hicParam.intraFlag);
    H = padconcatenation_sr(H,tempH,3);
end

%% perform matrix comparisons
% create data
rng(1)
N = 40;
T = 2;
A = zeros(N, N, T);
for i = 1:T
    temp = rand(N);
    A(:,:,i) = triu(temp)+triu(temp,1)';
end

%% Hi-C spector
% "The bin size was set to be 40kb"
% "The intra-chromosomal contact maps were not normalized (balanced)"
[Sd,Q] = hicSpector(A,B,20);

%% EPCA
% "17 Mb section of chr10, we binned the contacts into 250 kb bins"
% "We normalized these matrices by the expected number of random contacts at each genomic distance"
[epcaOut] = epca(A);

%% HiCRep
% smooth data
AHicRep = A;
hmat = [20  10;...
        11  25;...
        5   40;...
        3   100;...
        1   500;...
        0   1000];
h = 1;                  % smoothing filter span size h
hSmoothMat = ones(2*h+1)/(2*h+1)^2;
AHicRepSmooth = imfilter(AHicRep,hSmoothMat,'replicate');

% number check
figure, subplot(1,2,1), imagesc(A(:,:,1)), subplot(1,2,2), imagesc(AHicRepSmooth(:,:,1))
linkaxes(get(gcf,'children'))
sum(sum(A(:,:,1))), sum(sum(AHicRepSmooth(:,:,1)))

% remove counts > 5Mb from diagonal
binSizeMb = 1;
remMat = triu(ones(N),round(5*(1/binSizeMb)));
remMat = remMat+remMat';
AHicRepSmooth(find(repmat(remMat,[1 1 T]))) = NaN;

% Stratum-adjusted correlation coefficient (SCC)
kTotal = sum(~isnan(AHicRepSmooth(1,:,1)));
r1k = zeros(kTotal,1);
r2k = zeros(kTotal,1);
pk = zeros(kTotal,1);

for iK = 1:kTotal
    
    % get stratum contacts
    Xk = diag(AHicRepSmooth(:,:,1),iK-1);
    Yk = diag(AHicRepSmooth(:,:,2),iK-1);
    Nk = length(Xk);
    
    % 
    r1k(iK) = mean(Xk.*Yk)-(mean(Xk)*mean(Yk));
    r2k(iK) = sqrt(var(Xk)*var(Yk));
    r2k(iK) = sqrt(var(Xk)*var(Yk));
    
    % stratum-specific correlation ?k
    pk(iK) = r1k(iK)/r2k(iK);
end

% stratum-adjusted correlation coefficient (SCC)
ps = [];

%% EXTRA

% "several different resolutions (i.e., 10Kb, 25Kb, 40Kb, 100Kb, 500Kb, 1Mb)"
% "For most analysis, we used 40kb bins"
% "we chose to apply our method directly to raw data without bias correction"
%
% "only used the contacts within the range of 0-5Mb in the reproducibility assessment"
% "interactions over 5Mb in distance are rare (< 5% of reads) and highly stochastic"
% "Our method is applicable to both raw and bias corrected data"

% remove counts > 5Mb from diagonal
AHicRep = A;
binSizeMb = 1;
remMat = triu(ones(N),round(5*(1/binSizeMb)));

AHicRep(find(repmat(remMat,[1 1 T]))) = NaN;



% "The first stage is smoothing the raw contact matrix"
%
% "applying a 2D mean filter, which replaces the read count of each
% contact in the contact map with the average counts of all contacts in its
% neighborhood"

% "the optimal h should be adaptively chosen from the data"
% "We obtained h=20, 11, 5, 3, 1, and 0 for the resolution of 10Kb, 25Kb,
% 40Kb, 100Kb, 500Kb and 1Mb, respectively"
hmat = [20  10;...
        11  25;...
        5   40;...
        3   100;...
        1   500;...
        0   1000];

h = 1;      % smoothing filter span size h
hSmoothMat = ones(2*h+1)/(2*h+1)^2;
AHicRepSmooth = imfilter(AHicRep,hSmoothMat);
[Y,W] = ndnanfilter(AHicRep,hSmoothMat);


% "second stage, we apply a stratification approach to account for the
% pronounced distance dependence in the Hi-C data"
% 
% First we stratify the smoothed chromatin interactions according to their
% genomic distance



% then we apply a novel stratum-adjusted correlation coefficient statistic
% (SCC) to assess the reproducibility of the Hi-C matrices








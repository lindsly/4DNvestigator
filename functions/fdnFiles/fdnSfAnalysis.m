function [features,score] = fdnSfAnalysis(dataInfo,R,H,sampleSelect,goi,hExtract)
%fdnSfAnalysis Performs the 4DN Feature Analyzer
%   The 4DN Feature Analyzer is a method to analyze Hi-C and RNA-seq data
%   simultaneously. Centrality features from Hi-C are concatenated with
%   RNA-seq to form a feature matrix, which is then projected to a low
%   dimensional space. Details are describes in "Genome Architecture
%   Mediates Transcriptional Control of Human Myogenic Reprogramming"
%   doi: https://doi.org/10.1016/j.isci.2018.08.002
%
%   Input
%   dataInfo:       Data structure with experiment sample info
%   H:              Data structure containing all Hi-C data
%   R:              Data structure containing all RNA-seq data
%   sampleSelect:   Samples selected for analysis
%   goi:            Genes sets to be analyzed
%   hExtract:       Hi-C contact extraction method; 1MB ('1mb'), 100kb
%                   ('100kb'), or gene-level ('gene').
%
%   Output
%   features:       Feature matrix (e.g centrality and RNA-seq matrix)
%   score:          4DN score; location of bins in low dimensional
%                   projection
%
%   Scott Ronquist, scotronq@umich.edu. 2/3/19

%% Set default parameters
% Hi-C contact are extracted at 1Mb resolution by default
if ~exist('hExtract','var') || isempty(hExtract); hExtract = '1mb'; end

%% Select sample(s)
% select which samples will be analyzed
if ~exist('sampleSelect','var') || isempty(sampleSelect)
    hicSamples = dataInfo.sampleInfo.uniqueName(ismember(dataInfo.sampleInfo.dataType,'hic'));
    [tpIdx,~] = listdlg('PromptString','Select Samples to analyze',...
        'ListString',hicSamples,'ListSize',[500 500]);
    sampleSelect = hicSamples(tpIdx);
end

%% select GOI list
if ~exist('goi','var') || isempty(goi)
    [hallmarkGsea] = readMSigDBHallmark;
    hallmarkGseaNames = fieldnames(hallmarkGsea);
    
    [indx,~] = listdlg('PromptString','Select Genes to analyze',...
        'ListString',hallmarkGseaNames,'ListSize',[500 500]);
    
    goi = {};
    for iSelect = 1:length(indx)
        goi = [goi;hallmarkGsea.(hallmarkGseaNames{indx(iSelect)})];
    end
    goi = intersect(R.TPM.geneName,goi);
end

%% extract data
goiH = zeros(length(goi),length(goi),length(tpIdx));
goiR = zeros(length(goi),length(tpIdx));

for iSample = 1:length(tpIdx)
    
    % get RNA-seq column loc
    SampleTpRLocs = contains(R.TPM.Properties.VariableNames,...
        sprintf('%s_T%i',dataInfo.sampleInfo.sample{tpIdx(iSample)},...
        dataInfo.sampleInfo.timePoint(tpIdx(iSample))));
    
    % get RNA-seq values, replicate mean of TPM
    goiR(:,iSample) = mean(R.TPM{ismember(R.TPM.geneName,goi),SampleTpRLocs},2);
    
    % get Hi-C roiLoc, relative for extraction method
    switch hExtract
        case '1mb'
            binSize = 1E6;
            goiLoc = R.TPM{ismember(R.TPM.geneName,goi),1:3};
            goiLoc(:,2:3) = ceil(goiLoc(:,2:3)./binSize);
            
        case 'gene'
            goiLoc = R.TPM{ismember(R.TPM.geneName,goi),1:3};
            flankSize = 1E5;
            goiLoc(:,2) = goiLoc(:,2)-flankSize;
            goiLoc(:,3) = goiLoc(:,3)+flankSize;
    end
    
    % get Hi-C
    for iGoi = 1:length(goi)
        for iiGoi = iGoi:length(goi)
            disp([iGoi iiGoi])
            
            switch hExtract
                case '1mb'
                    goiH(iGoi,iiGoi,iSample) = mean(mean(H.s1mb.oe(...
                        H.s1mb.chrStart(goiLoc(iGoi,1))-1+goiLoc(iGoi,2):...
                        H.s1mb.chrStart(goiLoc(iGoi,1))-1+goiLoc(iGoi,3),...
                        H.s1mb.chrStart(goiLoc(iiGoi,1))-1+goiLoc(iiGoi,2):...
                        H.s1mb.chrStart(goiLoc(iiGoi,1))-1+goiLoc(iiGoi,3),...
                        tpIdx(iSample))));
                case 'gene'
                    error('In Progress, currently takes too long')
                    
                    % find 4DNvestigator juicertools path
                    juicerJarDir = mfilename('fullpath');
                    
                    juicerJarDirLevels = strfind(juicerJarDir,filesep);
                    juicerJarDir = [juicerJarDir(1:juicerJarDirLevels(end-1)),...
                        fullfile('hic','juicer2matlab','juicer_tools.jar')];
                    
                    % get Hi-C info
                    fnOut = 'juicer_temp.txt';
                    [status,cmdout] =...
                        system(sprintf('java -jar "%s" dump observed kr %s %i:%i:%i %i:%i:%i BP %i',...
                        juicerJarDir, SampleTpHPath{iSample},...
                        goiLoc(iGoi,1), goiLoc(iGoi,2), goiLoc(iGoi,3),...
                        goiLoc(iiGoi,1), goiLoc(iiGoi,2), goiLoc(iiGoi,3),...
                        min(dataInfo.hicHeader.BasePairdelimitedResolutions),fnOut));
                    
                    % report an error if status~=0
                    if status~=0
                        error(cmdout)
                    end
                    
                    % read data
                    temp = readtable(fnOut,'Delimiter','\t');
                    delete(fnOut)
                    
                    % add data to gene-level Hi-C, if available
                    if ~isempty(temp)
                        goiH(iGoi,iiGoi) = sum(temp{:,3});
                    end
            end
        end
    end
    goiH(:,:,iSample) = goiH(:,:,iSample) + triu(goiH(:,:,iSample),1)';
    
end

%% run sfAnalysis
[features,score] = sfAnalysis(goiH,goiR,goi);

end


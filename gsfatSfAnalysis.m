function [features,score] = gsfatSfAnalysis(dataInfo,R,H,sampleSelect,goi,hExtract)
%gsfatSfAnalysis Summary of this function goes here
%   Detailed explanation goes here
%
%   Input
%   dataInfo: 
%   R: 
%   H: 
%   sampleSelect: 
%   goi: 
%   hExtract: 
%
%   Output
%   features: 
%   score: 
%
%   Scott Ronquist, scotronq@umich.edu. 12/20/18

%% set default parameters
if ~exist('hExtract','var') || isempty(hExtract); hExtract = '1mb'; end

%% select sample(s)
if ~exist('sampleSelect','var') || isempty(sampleSelect)
    samples = unique(dataInfo.sampleInfo.sample);
    [indx,~] = listdlg('PromptString','Select Genes to analyze',...
        'ListString',samples,...
        'ListSize',[500 500],...
        'SelectionMode','single');
    sampleSelect = samples{indx};
end

for iSample = 1:length(sampleSelect)
    tps = unique(dataInfo.sampleInfo.timePoint);
    SampleTpRLocs = cell(length(tps),1);
    SampleTpHLocs = zeros(length(tps),1);
    for iTp = 1:length(tps)
        SampleTpRLocs{iTp} = ismember(R.TPM.Properties.VariableNames,...
            dataInfo.sampleInfo.name(not(cellfun('isempty',...
            regexp(dataInfo.sampleInfo.uniqueName,...
            sprintf('rnaseq_s%s_t%i*',sampleSelect,tps(iTp)))))));
        
        SampleTpHLocs(iTp) = dataInfo.sampleInfo.index(...
            not(cellfun('isempty',...
            regexp(dataInfo.sampleInfo.uniqueName,...
            sprintf('hic_s%s_t%i*',sampleSelect,tps(iTp))))));
        
        SampleTpHPath{iTp} = dataInfo.sampleInfo.path{...
            not(cellfun('isempty',...
            regexp(dataInfo.sampleInfo.uniqueName,...
            sprintf('hic_s%s_t%i*',sampleSelect,tps(iTp)))))};
    end
    
    %% select GOI list
    if ~exist('goi','var') || isempty(goi)
        [hallmarkGsea] = readMSigDBHallmark;
        hallmarkGseaNames = fieldnames(hallmarkGsea);
        
        [indx,tf] = listdlg('PromptString','Select Genes to analyze',...
            'ListString',hallmarkGseaNames,...
            'ListSize',[500 500]);
        
        goi = {};
        for iSelect = 1:length(indx)
            goi = [goi;hallmarkGsea.(hallmarkGseaNames{indx(iSelect)})];
        end
        goi = unique(goi);
    end
    
    %% extract goiH and goiR data
    % get gene locations, relative to binSize Hi-C
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
    
    % get Hi-C from gene list
    goiH = zeros(length(goi),length(goi),length(tps));
    goiR = zeros(length(goi),length(tps));
    for iTp = 1:length(tps)
        for i = 1:length(goi)
            goiR(i,iTp) = mean(R.TPM{ismember(R.TPM.geneName,goi{i}),SampleTpRLocs{iTp}},2);
            for ii = i:length(goi)
                disp([i ii])
                
                switch hExtract
                    case '1mb'
                        goiH(i,ii,iTp) = mean(mean(H.s1mb.oe(...
                            H.s1mb.chrStart(goiLoc(i,1))-1+goiLoc(i,2):...
                            H.s1mb.chrStart(goiLoc(i,1))-1+goiLoc(i,3),...
                            H.s1mb.chrStart(goiLoc(ii,1))-1+goiLoc(ii,2):...
                            H.s1mb.chrStart(goiLoc(ii,1))-1+goiLoc(ii,3),...
                            SampleTpHLocs(iTp))));
                    case 'gene'
                        % get Hi-C info
                        fnOut = 'juicer_temp.txt';
                        [status,cmdout] = system(sprintf('java -jar "%s" dump observed kr %s %i:%i:%i %i:%i:%i BP %i',...
                            'E:\MATLAB\gsfat\functions\juicer2matlab\juicer_tools.jar',...
                            SampleTpHPath{iTp},...
                            goiLoc(i,1),goiLoc(i,2),goiLoc(i,3),...
                            goiLoc(ii,1),goiLoc(ii,2),goiLoc(ii,3),...
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
                            goiH(i,ii) = sum(temp{:,3});
                        end
                end
                
            end
        end
        goiH(:,:,iTp) = max(cat(3,goiH(:,:,iTp),goiH(:,:,iTp)'),[],3);
    end
    
    %% run sfAnalysis
    [features{iSample},score{iSample}] = sfAnalysis(goiH,goiR,goi);
end

end


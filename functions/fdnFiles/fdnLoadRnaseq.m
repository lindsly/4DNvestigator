function [R] = fdnLoadRnaseq(dataInfo,H)
%fdnLoadRnaseq Loads and formats RNA-seq data specified in dataInfo
%   
%   Inputs
%   dataInfo:   structure containing all experiment metadata
%   H:          Hi-C data structure
%
%   Outputs
%   R:          RNA-seq data structure
%
%   Scott Ronquist, scotronq@umich.edu. 12/19/18

%% get information dataInfo
% get chr information from hic header
chrInfo = dataInfo.hicHeader.Chromosomes;
numChr = height(chrInfo);

% get Hi-C sample Locations
repGrp = findgroups(strcat(dataInfo.sampleInfo.sample(strcmp(dataInfo.sampleInfo.dataType,'rnaseq')),...
    num2str(dataInfo.sampleInfo.timePoint(strcmp(dataInfo.sampleInfo.dataType,'rnaseq')))));

% RSEM to mat
R = rsem2mat(dataInfo.sampleInfo.path(ismember(dataInfo.sampleInfo.dataType,'rnaseq')),...
    dataInfo.sampleInfo.name(ismember(dataInfo.sampleInfo.dataType,'rnaseq')),...
    dataInfo.hicHeader.refGenome);

Rfields = fields(R);

for i = 1:length(Rfields)
    R.(Rfields{i})(~cellfun(@isnumeric,R.(Rfields{i}).chr),:) = [];
    R.(Rfields{i}).chr = cell2mat(R.(Rfields{i}).chr);
end

%% bined RNA-seq
%100kb bins
chrBinSizes = ceil(chrInfo.chrLength/1E5);
[R.s100kb.tpm,R.s100kb.gene] = rna2bin(R.TPM{:,7:end},R.TPM.geneName,...
    [R.TPM.chr R.TPM.geneStart R.TPM.geneEnd],1E5,chrBinSizes);

% get mean for bins
R.s100kb.tpmMean = cell(numChr,1);
for iChr = 1:numChr
    R.s100kb.tpmMean{iChr} = zeros(size(R.s100kb.tpm{iChr},1),max(repGrp));
    for iRep = 1:max(repGrp)
        R.s100kb.tpmMean{iChr}(:,iRep) = mean(R.s100kb.tpm{iChr}(:,repGrp==iRep),2);
    end
end

%1mb bins
chrBinSizes = ceil(chrInfo.chrLength/1E6);
[R.s1mb.tpm,R.s1mb.gene] = rna2bin(R.TPM{:,7:end},R.TPM.geneName,...
    [R.TPM.chr R.TPM.geneStart R.TPM.geneEnd],1E6,chrBinSizes);

R.s1mb.tpm = vertcat(R.s1mb.tpm{:});
R.s1mb.gene = vertcat(R.s1mb.gene{:});

% get mean for bins
R.s1mb.tpmMean = zeros(size(R.s1mb.tpm,1),max(repGrp));
for iRep = 1:max(repGrp)
    R.s1mb.tpmMean(:,iRep) = mean(R.s1mb.tpm(:,repGrp==iRep),2);
end

%trim based on Hi-C
for iChr = 1:numChr
    R.s100kb.tpmTrim{iChr,1} = R.s100kb.tpm{iChr}(~H.s100kb.oeTrimBadLocs{iChr},:);
    R.s100kb.tpmMeanTrim{iChr,1} = R.s100kb.tpmMean{iChr}(~H.s100kb.oeTrimBadLocs{iChr},:);
    R.s100kb.geneTrim{iChr,1} = R.s100kb.gene{iChr}(~H.s100kb.oeTrimBadLocs{iChr},:);
end

R.s1mb.tpmTrim = R.s1mb.tpm(~H.s1mb.oeTrimBadLocs,:);
R.s1mb.tpmMeanTrim = R.s1mb.tpmMean(~H.s1mb.oeTrimBadLocs,:);
R.s1mb.geneTrim = R.s1mb.gene(~H.s1mb.oeTrimBadLocs,:);

end


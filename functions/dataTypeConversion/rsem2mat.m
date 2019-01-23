function [dataOut] = rsem2mat(fn,dataName,refGenome)
%rsem2mat takes the rsem raw data and formats it for MATLAB
%
%   Input
%   fn:         cell array of file names for RNA-seq samples
%   dataName:   names for each sample
%   refGenome:  references genome information
%
%   Output
%   dataOut:    formated RNA-seq data
%
%   Scott Ronquist, scotronq@umich.edu. 1/22/19

%% set default parameters
if ~exist('refGenome','var') || isempty(refGenome); refGenome='hg19'; end

% use file names if no specified names given
if ~exist('dataName','var') || isempty(dataName)
    dataName = cell(length(fn),1);
    for iSample = 1:length(fn)
        slashLoc = regexp(fn{iSample},'[\\/]');
        temp = fn{iSample}(slashLoc(end)+1:end);
        dotLoc = regexp(temp,'\.');
        dataName{iSample} = temp(1:dotLoc(1)-1);
    end
end

%% load geneinfo
switch refGenome
    case 'hg19'
        biomart = readtable('mart_export_ensembl_hg37_info.txt');
    case 'hg38'
        biomart = readtable('mart_export_ensembl_hg38_info.txt');
end

%% load and compare biomart info
for iSample = 1:length(fn)
    fprintf('formatting sample %d of %d\n',iSample,length(fn))
    
    dataIn = readtable(fn{iSample},'filetype', 'text');
    
    if iSample==1
        % keep genes that are found in both RNA-seq and biomart
        [C,~,IB] = intersect(dataIn.gene_id,biomart.GeneStableID);
        if ~isequal(dataIn.gene_id,C)
            error('biomart-RSEM gene incompatible')
        end
        geneName = biomart.HGNCSymbol(IB);
        geneDescription = biomart.GeneDescription(IB);
        chr = biomart.Chromosome_scaffoldName(IB);
        geneStart = biomart.GeneStart_bp_(IB);
        geneEnd = biomart.GeneEnd_bp_(IB);
        geneStrand = biomart.Strand(IB);
        
        % create table of gene information
        baseTable = [table(chr) table(geneStart) table(geneEnd),...
            table(geneName) table(geneDescription) table(geneStrand)]; % like .bed file (sort-of)
        
        dataOut.expected_count = baseTable;
        dataOut.TPM = baseTable;
        dataOut.FPKM = baseTable;
    end
    
    if ~isequal(dataIn.gene_id,C)
        error('biomart-RSEM gene incompatible')
    end
    
    dataOutFields = fields(dataOut);
    
    % add to tables with RNA-seq sample information
    for ii = 1:length(dataOutFields)
        dataOut.(dataOutFields{ii}) = [dataOut.(dataOutFields{ii}),table(dataIn.(dataOutFields{ii}))];
        dataOut.(dataOutFields{ii}).Properties.VariableNames{end} = dataName{iSample};
    end
    
end

for iField = 1:length(dataOutFields)
    fprintf('formatting RNAseq field: %s\n',dataOutFields{iField})
    
    %remove empty gene names
    dataOut.(dataOutFields{iField})(cellfun(@isempty,dataOut.(dataOutFields{iField}).geneName),:) = [];
    
    % chr str2num
    chrNums = cellfun(@str2num,dataOut.(dataOutFields{iField}).chr,'UniformOutput',0);
    dataOut.(dataOutFields{iField}).chr(~cellfun(@isempty,chrNums)) = chrNums(~cellfun(@isempty,chrNums));
    
    % sort
    dataOut.(dataOutFields{iField}) = sortrows(dataOut.(dataOutFields{iField}),'geneName');
end

end


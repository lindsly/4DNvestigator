function [dataOut] = rsem2mat(dataIn,refGenome)
%rsem2mat takes the rsem raw data and formats it for MATLAB
%   dataIn: table or filename or directory with RSEM output
%   dataOut: table for RSEM output, formatted with biomart

if nargin<2; refGenome='hg19'; end

%% load data
switch class(dataIn)
    %formatted table input
    case 'table'
        if ismember('dataType',dataIn.Properties.VariableNames)
            dataIn(~ismember(dataIn.dataType,'rnaseq'),:) = [];
        end
        
        dataFn = cell(height(dataIn),1);
        for iSample = 1:height(dataIn)
            dataFn{iSample} = [dataIn.folder{iSample},'\',dataIn.name{iSample}];
        end
        dataName = dataIn.uniqname;
        
        %file location input
    case 'char'
        if exist(dataIn)==7
            dataLoc = rdir([dataIn,'/**/*.genes.results'],'',1);
            
            dataFn = cell(height(dataIn),1);
            dataName = cell(height(dataIn),1);
            for iSample = 1:height(dataIn)
                dataFn{iSample} = [dataLoc(iSample).folder,'/',dataLoc(iSample).name];
                dotLoc = strfind(dataLoc(iSample).name,'.');
                dataName{iSample} = dataLoc(iSample).name(1:dotLoc(1)-1);
                dataName{iSample} = strrep(dataName,'-','_');
            end
        else
            dataFn = {dataIn};
            dataName = {'a'};
            
        end
    case 'cell' % FIX THIS LATER
        dataFn = dataIn;
        dataName = {'CT_r1','CT_r2','CT_r3','KT_r1','KT_r2','KT_r3'};
end

%% load geneinfo
switch refGenome
    case 'hg19'
        biomart = readtable('mart_export_ensembl_hg37_info.txt');
    case 'hg38'
        biomart = readtable('mart_export_ensembl_hg38_info.txt');
end

%% load and compare biomart info
for iSample = 1:length(dataFn)
    fprintf('formatting sample %d of %d\n',iSample,length(dataFn))
    
    dataIn = readtable(dataFn{iSample},'filetype', 'text');
    
    if iSample==1
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
        
        %         baseTable = [table(geneName) table(geneDescription) table(chr),...
        %             table(geneStart) table(geneEnd) table(geneStrand)];
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
    
    for ii = 1:length(dataOutFields)
        dataOut.(dataOutFields{ii}) = [dataOut.(dataOutFields{ii}),table(dataIn.(dataOutFields{ii}))];
        dataOut.(dataOutFields{ii}).Properties.VariableNames{end}=dataName{iSample};
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

%% EXTRA
% [LIA,LOCB] = ismember(biomart.GeneStableID,data_in.gene_id);
% biomart = [biomart table(LOCB)];
% biomart(~LIA,:) = [];

% temp_start = grpstats(biomart,'LOCB','min','DataVars','GeneStart_bp_');
% temp_start = grpstats(biomart,'LOCB','max','DataVars','GeneEnd_bp_');

% if length(unique(data_in.gene_id)) ~= length(data_in.gene_id)
%     error('gene list not unique')
% end
%
% [C,~,IB] = intersect(data_in.gene_id,biomart.GeneStableID);
% if ~isequal(data_in.gene_id,C)
%     error('biomart-RSEM gene incompatible')
% end
%
% gene_name = biomart.HGNCSymbol(IB);
% gene_description = biomart.GeneDescription(IB);
% chr = biomart.Chromosome_scaffoldName(IB);
% gene_start = biomart.GeneStart_bp_(IB);
% gene_end = biomart.GeneEnd_bp_(IB);
% gene_strand = biomart.Strand(IB);
%
% data_out = [table(gene_name) table(gene_description) table(chr),...
%     table(gene_start) table(gene_end) table(gene_strand) data_in];
%
% %remove empty gene names
% data_out(cellfun(@isempty,data_out.gene_name),:) = [];
%
% % chr str2num
% chr_nums = cellfun(@str2num,data_out.chr,'UniformOutput',0);
% data_out.chr(~cellfun(@isempty,chr_nums)) = chr_nums(~cellfun(@isempty,chr_nums));
%
% data_out = sortrows(data_out,'gene_name');
%
% %remove empty gene names
% data_out(cellfun(@isempty,data_out.gene_name),:) = [];
%
% % chr str2num
% chr_nums = cellfun(@str2num,data_out.chr,'UniformOutput',0);
% data_out.chr(~cellfun(@isempty,chr_nums)) = chr_nums(~cellfun(@isempty,chr_nums));
%
% % sort
% data_out = sortrows(data_out,'gene_name');

end


function [] = featureAnalyzerExample(Data_Loc, Folder_Result, chrSelect, dimReduc, binSize, sfType,featureType)
    %% 4DN Feature Analyzer example
    % This example shows how the "4DN feature analyzer" can be used to find
    % genes which change significantly in both structure and function over time
    %
    %   link to paper: (In preparation)
    %
    %   Version 1.0 (5/24/19)
    %   Written by: Scott Ronquist
    %   Contact:    scotronq@umich.edu
    %   Created:    5/24/19
    %
    %   Revision History:
    %   v1.0 (5/24/19)
    %   * featureAnalyzerExample.m created
    %
    %   v2.0 (10/21/20)
    %   * conversion to function, called from ExampleScript.m

    %% Load Data
    load(Data_Loc);
    load('ncbi_gene_info.mat');

    %% Set default 4DN Feature Analyzer parameters
    if ~exist('chrSelect','var')||isempty(chrSelect);chrSelect = 11;end
    if ~exist('dimReduc','var')||isempty(dimReduc);dimReduc='pca';end
    if ~exist('binSize','var')||isempty(binSize);binSize=1E5;end
    topEllipseFrac = .1;
    if ~exist('sfType','var')||isempty(sfType);sfType='sfmatrix';end
    if ~exist('featureType','var')||isempty(sfType);sfType='rna';end

    
    %% Extract region to analyze from Hi-C and RNA-seq
    switch binSize
        case 1E6
            goiH = H.s1mb.oeTrim(H.s1mb.chrStartTrim(chrSelect):H.s1mb.chrStartTrim(chrSelect+1)-1,...
                H.s1mb.chrStartTrim(chrSelect):H.s1mb.chrStartTrim(chrSelect+1)-1,:);
            if strcmp(featureType,'rna')
                goiR = R.s1mb.tpmMeanTrim(H.s1mb.chrStartTrim(chrSelect):H.s1mb.chrStartTrim(chrSelect+1)-1,:);
                goi = R.s1mb.geneTrim(H.s1mb.chrStartTrim(chrSelect):H.s1mb.chrStartTrim(chrSelect+1)-1);
            elseif strcmp(featureType,'other')
                % V Once other data is loaded, comment out the error below V
                error('Please load your data signal and gene information for 1 Mb bins');
                % Assign genes to goi in 1 Mb bins (rows) and time points (columns)
%                 goi = input_data_genes;
                % Assign 1D signal to goiR in 1 Mb bins
%                 goiR = input_data_signal;
            end
        case 1E5
            goiH = H.s100kb.oeTrim{chrSelect};
            if strcmp(featureType,'rna')
                goiR = R.s100kb.tpmMeanTrim{chrSelect};
                goi = R.s100kb.geneTrim{chrSelect};
            elseif strcmp(featureType,'other')
                error('Please load your data signal and gene information for 1 Mb bins');
                % Assign 1D signal to goiR in 100 kb bins (rows) and time points (columns)
%                 goiR = input_data_signal;
                % Assign genes to goi in 100 kb bins
%                 goi = input_data_genes;
            end
        otherwise
            error('please select a correct bin size (1E6 or 1E5)')
    end

    %% Run the 4DNfeature analyzer visualization
    [features,score,genes] = sfAnalysis(goiH,goiR,goi,[],[],dimReduc,topEllipseFrac,sfType,featureType);
    
    %% Format genes in loci with largest structure-function changes
    genes = unique(genes,'stable');
    genes = intersect(genes, ncbi_gene_info.Symbol,'stable');
    ID=ncbi_gene_info.GeneID(genes);
    for i = 1:length(genes)
        url_ncbi{i,1}=sprintf("https://www.ncbi.nlm.nih.gov/gene/%d",ID(i));
        url_genecards{i,1}=sprintf("https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s",genes{i});
    end

    idx = ismember(R.TPM.geneName,genes);
    genes_table = R.TPM(idx,1:5);
    [~,order] = ismember(genes, genes_table.geneName);
    genes_table = genes_table(order,:);
    genes_table = movevars(genes_table,'geneName','Before','chr');
    genes_table.Properties.RowNames = genes_table.geneName;
    genes_table.NCBI_web = url_ncbi;
    genes_table.GeneCards_web = url_genecards;

    %% Save figures
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
      FigHandle = FigList(iFig);
      FigName   = get(FigHandle, 'Name');
      savefig(FigHandle, [Folder_Result, '\',FigName, '.fig']);
      saveas(FigHandle, [Folder_Result, '\',FigName, '.png']);
    end
end
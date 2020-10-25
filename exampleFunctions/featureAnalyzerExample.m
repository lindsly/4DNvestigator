function [] = featureAnalyzerExample(Data_Loc, Folder_Result, chrSelect, dimReduc, binSize)
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
    % 4DNvestigator formatted Hi-C and RNA-seq data is available here:
    %   https://drive.google.com/open?id=11XQ0CudxRvM5P6aI57yes2h8XXLhxlRr

    % Raw data (.hic and .results) files can be found here:
    %   https://drive.google.com/open?id=1lSyU-7I0ME3X70Mt_-HjLMtPc-BMKHxm

    % Add 4DNvestigator tools to path
%     filepath = mfilename('fullpath');
%     fdnPath = filepath(1:strfind(filepath,'4DNvestigator')+12);
%     addpath(genpath(fdnPath))

    % Load data
    load(Data_Loc);

    %% 4DN Feature Analyzer Example
    % Select 4DN Feature Analyzer parameters
    % This example is set up to analyze an entire genome at 1 Mb or 100 kb
    % resolution

    %% Set default 4DN Feature Analyzer parameters
    if ~exist('chrSelect','var')||isempty(chrSelect);chrSelect = 11;end
    if ~exist('dimReduc','var')||isempty(dimReduc);dimReduc='pca';end
    if ~exist('binSize','var')||isempty(binSize);binSize=1E5;end
    topEllipseFrac = .1;
    
    %% Extract region to analyze from Hi-C and RNA-seq
    switch binSize
        case 1E6
            goiH = H.s1mb.oeTrim(H.s1mb.chrStartTrim(chrSelect):H.s1mb.chrStartTrim(chrSelect+1)-1,...
                H.s1mb.chrStartTrim(chrSelect):H.s1mb.chrStartTrim(chrSelect+1)-1,:);
            goiR = R.s1mb.tpmMeanTrim(H.s1mb.chrStartTrim(chrSelect):H.s1mb.chrStartTrim(chrSelect+1)-1,:);
            goi = R.s1mb.geneTrim(H.s1mb.chrStartTrim(chrSelect):H.s1mb.chrStartTrim(chrSelect+1)-1);
        case 1E5
            goiH = H.s100kb.oeTrim{chrSelect};
            goiR = R.s100kb.tpmMeanTrim{chrSelect};
            goi = R.s100kb.geneTrim{chrSelect};
        otherwise
            error('please select a correct bin size (1E6 or 1E5)')
    end

    %% Run the 4DNfeature analyzer
    [features,score] = sfAnalysis(goiH,log2(goiR+1),goi,[],[],[],dimReduc,topEllipseFrac);

    %% Save figure
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
      FigHandle = FigList(iFig);
      FigName   = get(FigHandle, 'Name');
      savefig(FigHandle, [Folder_Result, '\',FigName, '.fig']);
    end
    
end
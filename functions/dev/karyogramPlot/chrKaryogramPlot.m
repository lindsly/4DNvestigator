function [h] = chrKaryogramPlot(geneList,rnaSeq)
%chrKaryogramPlot Plots the human genome karyogram with specified gene
%positions label. Color represents differential expression
%
%   Input:
%   geneList: list of genes to be labeled
%   rnaSeq: expression values, typically log2FC
%
%   Output:
%   h: figure handle
% 
%   Scott Ronquist, 11/16/18

%% Plot karyogram
hs_cytobands = cytobandread('hs_cytoBand.txt');
[h,appdata] = chromosomeplotSR(hs_cytobands); hold on

% get chr plot info from appdata
chrLengthBp = [findpeaks(double(hs_cytobands.BandEndBPs));double(hs_cytobands.BandEndBPs(end))];
chrLengthPlot = appdata.chr_len;
chrXLim = appdata.chr_xlim;
chrYLim = appdata.chr_ylim;
chrLabels = appdata.chromLabels;

% load gene position information
bioMart = readtable('mart_export_ensembl_hg37_info.txt');
[~,IA,~] = unique(bioMart.HGNCSymbol);
bioMart = bioMart(IA,:);

%% get color scale info from rnaSeq
% create a default color map ranging from red to green
cLength = 100;
red = [1 0 0];
green = [0 1 0];
colorScale = [linspace(red(1),green(1),cLength)',...
    linspace(red(2),green(2),cLength)',...
    linspace(red(3),green(3),cLength)'];
rnaSeqNorm = ceil(((rnaSeq - min(rnaSeq))/(max(rnaSeq) - min(rnaSeq)))*cLength); % normalize 0 to cLength
rnaSeqNorm(rnaSeqNorm==0) = 1;

%% Plot genes
missingGene = {'MISSING GENES';''};
geneBuff = 1E6;     % add length to gene to visualize better
for iGene = 1:length(geneList)
    % get gene loc
    tempGeneBioMart = bioMart(ismember(bioMart.HGNCSymbol,geneList{iGene}),:);
    
    % if gene not found, list here
    if isempty(tempGeneBioMart)
        missingGene = [missingGene;geneList{iGene}];
        continue
    end
    
    % convert chr Bp coords to Figure coords
    tempChr = find(ismember(chrLabels,tempGeneBioMart.Chromosome_scaffoldName));
    
    % if gene not found, list here
    if isempty(tempChr)
        missingGene = [missingGene;geneList{iGene}];
        continue
    end
    
    tempPlotGeneLoc = chrLengthPlot(tempChr)-...
        (([tempGeneBioMart.GeneStart_bp_-geneBuff,tempGeneBioMart.GeneEnd_bp_+geneBuff]/chrLengthBp(tempChr))*...
        chrLengthPlot(tempChr));
    
    % Plot gene position
    tempX = chrXLim(:,tempChr);
    tempY = chrYLim(:,tempChr);
    tempColor = colorScale(rnaSeqNorm(iGene),:);
    
    patch([tempX(1) tempX(1) tempX(2) tempX(2)],...
        [tempY(1)+tempPlotGeneLoc(1), tempY(1)+tempPlotGeneLoc(2),...
        tempY(1)+tempPlotGeneLoc(2), tempY(1)+tempPlotGeneLoc(1)],...
        tempColor, 'edgecolor', 'none')
    text(tempX(2), tempY(1)+tempPlotGeneLoc(1), geneList{iGene})
end

%% Plot genes not found in top right
text(.65, .95, missingGene)

%% Plot colorbar
axes('position',[.9 .85 .03 .12]), image(reshape(colorScale,[size(colorScale,1),1,3]))
set(gca,'xtick',[],'ytick',[1 size(colorScale,1)],...
    'YAxisLocation', 'right','yticklabels',[round(min(rnaSeq),2) round(max(rnaSeq),2)])

end


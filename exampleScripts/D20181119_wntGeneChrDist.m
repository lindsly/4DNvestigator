% Goal: determine the distribution of WNT genes on chromosomes
clear
close all

%% karyogram with GOIs and expression
load('GeneTADinfo.mat', 'GeneInfo')
geneLoc = cell2mat(GeneInfo(:,4));

%load WNT genes
wntKeggGenes = readtable('wntKeggGenes.xlsx','ReadVariableNames',0);
wntGenesTargetMarkus = readtable('Wnt_Genes.xlsx','Sheet',3);
wntGenesAll = unique([wntKeggGenes{:,2};wntGenesTargetMarkus{:,1};wntGenesTargetMarkus{:,2}]);
wntGenesAll = wntGenesAll(2:end);% remove empty top row

chrKaryogramPlot(wntGenesAll,rand(length(wntGenesAll),1))

%% EXTRA load GeneInfo
% load('GeneTADinfo.mat', 'GeneInfo')
% geneLoc = cell2mat(GeneInfo(:,4));
% 
% %load WNT genes
% wntKeggGenes = readtable('wntKeggGenes.xlsx','ReadVariableNames',0);
% wntGenesTargetMarkus = readtable('Wnt_Genes.xlsx','Sheet',3);
% wntGenesAll = unique([wntKeggGenes{:,2};wntGenesTargetMarkus{:,1};wntGenesTargetMarkus{:,2}]);
% wntGenesAll = wntGenesAll(2:end);% remove empty top row
% 
% % get WNT gene locs
% wntGeneLoc = ismember(GeneInfo(:,1),wntGenesAll);
% 
% %% get percentage of WNT genes
% wntChrPrct = [];
% for iChr = 1:22
%     wntChrPrct(iChr) = (sum(geneLoc(wntGeneLoc,1) == iChr)/sum(geneLoc(:,1) == iChr))*100;
% end
% 
% % figure
% figure, bar(wntChrPrct)
% title('percentage of genes involved in WNT, per chr')
% 
% %% plot chromosome
% hs_cytobands = cytobandread('hs_cytoBand.txt');
% [h,appdata] = chromosomeplotSR(hs_cytobands); hold on
% chrLength = findpeaks(double(hs_cytobands.BandEndBPs));
% 
% for iWntGene = 1:length(wntGenesAll)
%     % get gene loc
%     tempGeneLoc = geneLoc(ismember(GeneInfo(:,1),wntGenesAll{iWntGene}),:);
%     if isempty(tempGeneLoc)
%         continue
%     end
%         
%     tempGeneLoc(2:3) = [tempGeneLoc(2)-1E6,tempGeneLoc(3)+1E6];
%     
%     % get chr coords
%     tempPlotGeneLoc = appdata.chr_len(tempGeneLoc(1))-...
%         ((tempGeneLoc(2:3)./chrLength(tempGeneLoc(1)))*appdata.chr_len(tempGeneLoc(1)));
%     tempX = appdata.chr_xlim(:,tempGeneLoc(1));
%     tempY = appdata.chr_ylim(:,tempGeneLoc(1));
%     
%     % add to figure
%     xBuff = .00;
%     patch([tempX(1)-xBuff tempX(1)-xBuff tempX(2)+xBuff tempX(2)+xBuff],...
%         [tempY(1)+tempPlotGeneLoc(1), tempY(1)+tempPlotGeneLoc(2),...
%         tempY(1)+tempPlotGeneLoc(2), tempY(1)+tempPlotGeneLoc(1)],...
%         [0 0 1], 'edgecolor', 'none')
%     text(tempX(2)+xBuff, tempY(1)+tempPlotGeneLoc(1), wntGenesAll{iWntGene})
% end
% 
% title('Human Karyogram')


%% EXTRA
% load coriell_baccgh
% Struct = cghfreqplot(coriell_data);
% S = cghcbs(coriell_data,'sampleindex',3,'showplot',true);
% 
% cnvStruct = struct('Chromosome', [10 11],...
%  'CNVType', [2 1],...
%  'Start', [S.SegmentData(10).Start(2),...
%   S.SegmentData(11).Start(2)]*1000,...
%  'End',   [S.SegmentData(10).End(2),...
%   S.SegmentData(11).End(2)]*1000)
% 
% chromosomeplot('hs_cytoBand.txt', 'CNV', cnvStruct, 'unit', 2)



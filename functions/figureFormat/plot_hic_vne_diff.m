function [h] = plot_hic_vne_diff(vne,Sample_names)
%plot_hic_vne_diff plots the difference in Hi-C VNE per chromosome
%   vne: VNE for each chromosome, 2 x 22 matrix
%   Sample_names: 2 x 1 cell with the names of sample for title
%
%   Scott Ronquist,6/15/18

h = figure('position',[100 100 800 400]);
b = bar(vne(1,:)-vne(2,:));
title(sprintf('%s - %s',strrep(Sample_names{1},'_',' '),strrep(Sample_names{2},'_',' ')))

set(gca,'LineWidth',2,'FontSize',15,'xtick',1:22,'xticklabels',1:22)
ylabel('VNE Difference')
xlabel('Chromosome')
axis tight

end


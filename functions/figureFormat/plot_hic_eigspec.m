function [h] = plot_hic_eigspec(eig_spec,Sample_name)
%plot_hic_eigspec plots the eigenspectrum derived from Hi-C
%   eig_spec: eigenspectrum for each chromosome, 22 x N matrix
%   Sample_name: name of sample for title
%
%   Scott Ronquist,6/15/18

temp_color = hsv(size(eig_spec,2));

h = figure('position',[100 100 800 400]);
b = bar(eig_spec,'stacked','FaceColor','flat');
for eig_color = 1:size(eig_spec,2)
    b(eig_color).CData = repmat(temp_color(eig_color,:),size(b(eig_color).CData,1),1);
end
title(sprintf('%s Normalized Eigenvalues',strrep(Sample_name,'_',' ')))
set(gca,'LineWidth',2,'FontSize',15,'xtick',1:22,'xticklabels',1:22)
ylabel('Normalized Eigen Spectrum')
xlabel('Chromosome')
axis tight
end


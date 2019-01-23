function [outputArg1,outputArg2] = plot_R_FDV_hic(R,FDV,H,R_norm,H_norm)
%plot_R_FDV_hic plots RNA-seq, Fiedler vector and Hi-C in same plot,
%aligned
%   R: RNA-seq vector
%   FDV: fiedler vector
%   H: Hi-C matrix
%   R_norm: [0 1], log2 normalization
%   H_norm: [0 1], Hi-C correlation matrix

if nargin < 4;R_norm=0;end
if nargin < 5;H_norm=0;end

if R_norm
    R = log2(R+.5)+1;
end

if H_norm
    H = corr(H);
    H = H-diag(diag(H));
end

figure('position',[100 100 700 1000]),
subplot(6,1,1)
bar(R)
xlim([0 length(R)])
ax = gca;ax.LineWidth = 2;ax.FontSize=8;
ylabel('log2(TPM)','fontsize',12)
title('RNA-seq','fontsize',15,'fontweight','bold')

subplot(6,1,2)
temp = FDV;
bar(find(temp<0),temp(temp<0),'r'),hold on
bar(find(temp>=0),temp(temp>=0),'g')
xlim([0 length(temp)])
ax = gca;ax.LineWidth = 2;ax.FontSize=8;
ylabel('normalized [-1,1]','fontsize',12)
title('Chromatin Patterning','fontsize',15,'fontweight','bold')

subplot(6,1,3:6)
imagesc(H)
ax = gca;ax.LineWidth = 2;ax.FontSize=8;
title('Hi-C','fontsize',15,'fontweight','bold')

% ax = suptitle(sprintf('%s, chr%d, 100kb',...
%     strrep(samples.names{sample_select},'_',' '),chr_select));
% ax.FontWeight='bold';
end


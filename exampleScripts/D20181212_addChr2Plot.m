% Goal: create a function that adds a chromosome to the side of an axis

%% load paths/defaults
clear
close all
addpath(genpath('E:\MATLAB\gsfat'))

%% plot Hi-C
%%% PARAMETERS vvv
hicParam.binType = 'BP';
hicParam.binSize = 1E5;
hicParam.norm1d = 'KR';
hicParam.norm3d = 'oe';
hicParam.intraFlag = 1;
hicParam.chr = 1;
hicParam.fn = 'E:\MATLAB\gsfat\sampleData\hic\aldh_N.hic';
%%% PARAMETERS ^^^

H = hic2mat(hicParam.norm3d,hicParam.norm1d,hicParam.fn,...
    hicParam.chr,hicParam.chr,hicParam.binType,hicParam.binSize,hicParam.intraFlag);

figure, imagesc(log(H))
axis square

%% plot chromosome
% needs to plot next to current figure axis
% needs to be able to trim chromosome image

hs_cytobands = cytobandread('hs_cytoBand.txt');
chromosomeplotSR('hs_cytoBand.txt', hicParam.chr, 'addtoplot', gca,...
    'Unit', hicParam.binSize);

%% end script
if 1==1
    return
end

%% EXAMPLE: subplot overlap attempt
a = rand(5);
figure, imagesc(a), axis square
colorbarScaled(.3)
hax = axes('Position', [.35, .35, .3, .3]);

%% EXAMPLE: subplot overlap example
figure
y = zeros(4,15);
for k = 1:4
    y(k,:) = rand(1,15);
    subplot(2, 2, k)
    plot(y(k,:));
end
hax = axes('Position', [.35, .35, .3, .3]);
bar(hax,y,'EdgeColor','none')
set(hax,'XTick',[])

%% EXAMPLE: chr addtoplot
hs_cytobands = cytobandread('hs_cytoBand.txt');

load coriell_baccgh
S = cghcbs(coriell_data,'sampleindex',3,'chromosome',10,...
    'showplot',10);
set(gcf,'color','w'); % Set the background of the figure to white.
chromosomeplotSR('hs_cytoBand.txt', 10, 'addtoplot', gca,...
    'Unit', 2);
       
       
       
       
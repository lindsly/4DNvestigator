function [h] = plot_hic(a,c_map_select,c_bar_switch)
%plot_hic Plot Hi-C matrices using Rajapakse Lab standards
%   a is the Hi-C matrix
%   c_map_select is the colormap style
%   c_bar_switch is 0 or 1 to turn colorbar off or on
%   h is the axis handle

%   default color map is inverted 'hot', used by Dekker
%   https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0831-x
%   'erez' colormap is white to red, used in 
%   https://www.sciencedirect.com/science/article/pii/S0092867414014974#fig2
%   Scott Ronquist, 4/21/18

if nargin<3; c_bar_switch=0;end
if nargin<2; c_map_select='dekker';end

%% format (to remove NaNs)
temp = a;
temp(isinf(temp)) = 0;
temp2 = sum(temp) == 0;
temp(temp2,:) = NaN;
temp(:,temp2) = NaN;

%% plot
h = imagesc(temp);
axis equal
xlim([.5 size(a,2)+.5])
ylim([.5 size(a,1)+.5])

h.Parent.LineWidth=2;
h.Parent.FontSize=15;
set(h,'AlphaData',~isnan(temp))

%% colormap
switch c_map_select
    case 'dekker'
        temp = colormap('hot');
        colormap(flipud(temp));
    case 'erez'
        temp = [ones(64,1),[1:-1/63:0]',[1:-1/63:0]'];
        colormap(temp);
end

if c_bar_switch
%     colorbar('position',[.87 .73 .05 .2]);
    temp = h.Parent.Position;
    colorbar('position',[temp(1)+temp(3) temp(2) .05 temp(4)]);
end

end


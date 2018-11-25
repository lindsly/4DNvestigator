function [] = addROICircles(A,circColor,circRadius,circLineWidth)
%addROICircles adds circles to regions of interest
%
%   Input:
%   figAxis: axis handle
%   A: Logical matrix with ROIs
%   circColor: circle color
%
%   Scott Ronquist, 11/16/18. scotronq@umich.edu

%% Start
if nargin < 2;circColor = 'blue';end
if nargin < 3;circRadius = 1.5;end
if nargin < 4;circLineWidth = 2;end

stats = regionprops('table',A,'Centroid',...
    'MajorAxisLength','MinorAxisLength');
centers = stats.Centroid;
diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
radii = diameters/2;
viscircles(centers,radii*circRadius,'Color',circColor,'lineWidth',circLineWidth);
end


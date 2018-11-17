function [] = addROICircles(A,circColor)
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

stats = regionprops('table',A,'Centroid',...
    'MajorAxisLength','MinorAxisLength');
centers = stats.Centroid;
diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
radii = diameters/2;
viscircles(centers,radii*1.5,'Color',circColor);
end


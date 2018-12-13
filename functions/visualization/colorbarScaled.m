function [cb] = colorbarScaled(cBarScale)
%colorbarScaled adds a colorbar to the current axis, at specified length
%
%   Input
%   cBarLength: percentage of colorbar length (default=.5)
%
%   Output
%   h: colorbar MATLAB handle
%
%   Example
%   a = rand(5);
%   figure, imagesc(a), axis square
%   colorbarScaled(.3)
%
%   Scott Ronquist, scotronq@umich.edu
%% default parameters
if ~exist('cBarLength','var')||isempty(cBarScale);cBarScale=.5;end

%% add scaled colorbar
% get axis position before colorbar is added
ax = gca;
pos2 = plotboxpos(ax);

% add colorbar
cb = colorbar;

% move colorbar
cb.Position = [pos2(1)+pos2(3)+.01,...
    pos2(2)+pos2(4)*(1-cBarScale),...
    cb.Position(3), pos2(4)*cBarScale];

%% extra
% ax = gca;
% [ax.OuterPosition;
% ax.TightInset;
% ax.Position]

end


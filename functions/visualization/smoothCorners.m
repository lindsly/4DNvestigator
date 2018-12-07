function [finalPlot,linearX] = smoothCorners(y,x,corners,interpNum)
%smoothCorners takes x,y data, finds corners, then smooths corners
%   following methods outlined in:
%   https://stackoverflow.com/questions/43720500/smooth-out-a-single-sharp-corner-in-a-function-in-matlab
%
%   Scott Ronquist, 11/26/18. scotronq@umich.edu

%% default values
corners = [];       % fix this later is user wants to define where corners are
interpNum = .01;    % fix this later
x = 1:length(y);    % fix this later

%% linear interpolate data
% first linear interp
linearX = x(1):.01:x(end);
linearY = interp1(x,y,linearX);

%% find corners
D2 = diff(linearY, 2);
roundMe = find(abs(D2)>0.1*max(D2))+1; % offset by 1 because second derivative

if 1==0
    figure;
    subplot(3,1,1)
    plot(linearY); title 'shark plot'
    hold on;
    plot(roundMe, linearY(roundMe),'r*')
    legend('input','corners found')
end

%% smooth 1
% take N points on either side of the sharp corners:
N = 3;

% boxplot filtered version of the curve
boxFilt = ones(1, 2*N+1)/(2*N+1);
smoothShark1 = convn(linearY, boxFilt, 'same'); % box plot

% second filter - smoother
smoothShark2 = convn(smoothShark1, boxFilt, 'same');

if 1==0
    subplot(3,1,2)
    plot(linearY)
    hold on
    plot(smoothShark1);
    hold on
    plot(smoothShark2);
end

%% filter
% Now apply filtering only to points near the discontinuity
smoothMe = zeros(size(linearY));
smoothMe(roundMe)=1;
smoothMe = convn(smoothMe, boxFilt, 'same');
smoothMe(smoothMe>0)=1; % this finds N points on either side of the corner

finalPlot=linearY;
smoothIndx = find(smoothMe);
finalPlot(smoothIndx)=smoothShark2(smoothIndx);

if 1==0
    subplot(3,1,3)
    plot(linearY)
    hold on
    plot(finalPlot,'g')
    plot(smoothIndx, finalPlot(smoothIndx), 'r*')
end
end


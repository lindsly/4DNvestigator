function [dataMean,dataStdLim] = smoothPlotMeanStd(data,interpNum)
%smoothPlotMeanStd plots data mean with errorbounds, interpolated
%
%   Input
%   data: is input data
%
%   Scott Ronquist, 11/19/18. scotronq@umich.edu
%% set defaults
if nargin < 2;interpNum = .1;end

%% get data mean and std
dataMean = mean(data,1);
dataStdLim = [dataMean+nanstd(data,1);dataMean-nanstd(data,1)];

% smooth data
dataMeanSmooth = interp1(1:size(dataMean,2),...
    dataMean,1:interpNum:size(dataMean,2),'spline');
dataStdLimSmooth = [interp1(1:size(dataStdLim(1,:),2),...
    dataStdLim(1,:),1:interpNum:size(dataStdLim(1,:),2),'spline');...
    interp1(1:size(dataStdLim(2,:),2),dataStdLim(2,:),...
    1:interpNum:size(dataStdLim(2,:),2),'spline')];

%% plot
patch([1:interpNum:size(dataMean,2),size(dataMean,2):-interpNum:1],...
    [dataStdLimSmooth(1,:),fliplr(dataStdLimSmooth(2,:))],'b',...
    'facealpha',.2,'edgecolor','none')
plot(1:interpNum:size(dataMean,2),dataMeanSmooth','b')
plot(dataMean,'bs')
axis tight

end


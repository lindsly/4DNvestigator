function [dataMean,dataStdLim] = smoothPlotMeanStd(data)
%smoothPlotMeanStd plots data mean with errorbounds, interpolated
%
%   Input
%   data: is input data
%
%   Scott Ronquist, 11/19/18. scotronq@umich.edu

% get data mean and std
dataMean = mean(data);
dataStdLim = [dataMean+nanstd(data);dataMean-nanstd(data)];

% smooth data
dataMeanSmooth = interp1(1:size(dataMean,2),...
    dataMean,1:.1:size(dataMean,2),'spline');
dataStdLimSmooth = [interp1(1:size(dataStdLim(1,:),2),...
    dataStdLim(1,:),1:.1:size(dataStdLim(1,:),2),'spline');...
    interp1(1:size(dataStdLim(2,:),2),dataStdLim(2,:),...
    1:.1:size(dataStdLim(2,:),2),'spline')];

%% plot
patch([1:.1:size(dataMean,2),size(dataMean,2):-.1:1],...
    [dataStdLimSmooth(1,:),fliplr(dataStdLimSmooth(2,:))],'b',...
    'facealpha',.2,'edgecolor','none')
plot(1:.1:size(dataMean,2),dataMeanSmooth','b')
plot(dataMean,'bs')
axis tight

end


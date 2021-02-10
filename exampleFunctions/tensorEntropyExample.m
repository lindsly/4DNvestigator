function [] = tensorEntropyExample(Data_Loc, Folder_Result)
    %% 4DN Feature Analyzer example
    % This example shows how to compute tensor entropy using
    % hicTensorEntropy.m and offers a comparison to von Neumann Entropy
    % 
    %   Written by: Can Chen, Stephen Lindsly
    %   Contact:    lindsly@umich.edu
    %   Created:    1/20/21
    % 
    %   Revision History:
    %   v1.0 (1/20/21)
    %   Script written to calculate tensor entropy for comparison to von
    %   Neumann entropy by Can Chen
    %
    %   v1.1 (1/27/21)
    %   Adapted into an example function called by ExampleScript.m and 
    %   moved tensor entropy calculation to the function hicTensorEntropy.m
    %   by Stephen Lindsly
    
    %% Load Data
    % Add 4DNvestigator tools to path
    filepath = mfilename('fullpath');
    fdnPath = filepath(1:strfind(filepath,'4DNvestigator')+12);
    addpath(genpath(fdnPath))

    % Load example MYOD reprogramming and proliferation data
    load(Data_Loc)
    A1 = chr14_myod;
    A2 = chr14_fib;

    % load('myodEntropy.mat');

    k = 3;

    combData = cat(3, chr14_fib, chr14_myod);
    thresh = 0.95;
    
    %% Calculate Tensor Entropy of Time Series Hi-C
    tensorEntropy_A1 = hicTensorEntropy(A1,thresh);
    tensorEntropy_A2 = hicTensorEntropy(A2,thresh);

    %% Calculate von Neumann network entropy for comparison
    entropyMatrix = zeros(2, size(combData, 3));
    for i = 1:size(combData, 3)
        % Network Entropy
        adjMatrix = combData(:, :, i);
        adjMatrix = adjMatrix - diag(diag(adjMatrix));
        adjMatrix = adjMatrix>=0.95;
        lapMatrix = diag(sum(adjMatrix, 1))-adjMatrix;
        eigValues = eig(lapMatrix);
        eigValues = eigValues/sum(eigValues);
        eigValues = eigValues(eigValues>0);
        entropyMatrix(1, i) = -sum(eigValues.*log(eigValues));
    end
   
    %% Plotting
    % Get colors for plots
    color_opts = color_options;
    
    % Plot tensor entropy
    figure('Position', [391 453 1023 420])
    subplot(1,2,1)
        data = tensorEntropy_A1;
        dataMeanSmooth = interp1(1:size(data,2),data,1:.1:size(data,2),'makima');
        p(1) = plot(1:.1:size(data,2),dataMeanSmooth,'Color',color_opts(1,:),'LineWidth',1.5); % Reprogramming
        hold on
        plot(1:size(data,2),data,'.','Color',color_opts(1,:),'MarkerSize',20)
        data = tensorEntropy_A2;
        dataMeanSmooth = interp1(1:size(data,2),data,1:.1:size(data,2),'makima');
        p(2) =  plot(1:.1:size(data,2),dataMeanSmooth,'Color',color_opts(2,:),'LineWidth',1.5); % Proliferation
        plot(1:size(data,2),data,'.','Color',color_opts(2,:),'MarkerSize',20)

        ylim_temp = get(gca,'YLim');
        xlim_temp = get(gca,'XLim');
%         set(gca,'LineWidth',1.25,'FontSize',12,'YLim',[ylim_temp(1)*.95 ylim_temp(2)*1.05],'xLim',[0 xlim_temp(2)+1]);
        title('Tensor Entropy')
        legend(p,{'Reprogramming','Proliferation'},'Location','southwest')
        axis square

    % Plot von Neumann entropy for comparison
    subplot(1,2,2)    
        data = entropyMatrix(1,1:8);
        dataMeanSmooth = interp1(1:size(data,2),data,1:.1:size(data,2),'makima');
        plot(1:.1:size(data,2),dataMeanSmooth,'Color',color_opts(2,:),'LineWidth',1.5) % Proliferation
        hold on
        plot(1:size(data,2),data,'.','Color',color_opts(2,:),'MarkerSize',20)
        data = entropyMatrix(1,9:16);
        dataMeanSmooth = interp1(1:size(data,2),data,1:.1:size(data,2),'makima');
        plot(1:.1:size(data,2),dataMeanSmooth,'Color',color_opts(1,:),'LineWidth',1.5) % Reprogramming
        plot(1:size(data,2),data,'.','Color',color_opts(1,:),'MarkerSize',20)
        
        ylim_temp = get(gca,'YLim');
        xlim_temp = get(gca,'XLim');
%         set(gca,'LineWidth',1.25,'FontSize',12,'YLim',[ylim_temp(1)*.95 ylim_temp(2)*1.05],'xLim',[0 xlim_temp(2)+1]);
        title('Network Entropy')
%         legend({'Proliferation','Reprogramming'})
        axis square

end

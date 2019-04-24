function [h] = plot_fdv(fdv,multi_plot)
%plot_fdv plots fiedler vector
%   fdv is the fiedler vector, this can have multiple rows (samples)
%   h is bar graph image handle
%   Scott Ronquist, 4/21/18

if nargin < 2; multi_plot=size(fdv,1)-1;end

for i = 1:size(fdv,1)
    %% plot bar graph
    if multi_plot; subplot(size(fdv,1),1,i); end
    
    h{i} = bar(fdv(i,:));
    xlim([.5 size(fdv,2)+.5])
    h{i}.FaceColor = 'flat';
    
    %% format colors, show switch
    pos_idx = fdv(i,:)>0;
    neg_idx = fdv(i,:)<0;
    
    h{i}.CData(pos_idx,:) = repmat([0 1 0],sum(pos_idx),1);
    h{i}.CData(neg_idx,:) = repmat([1 0 0],sum(neg_idx),1);
    
    if i>1
        switch_idx = fdv(i,:).*fdv(i-1,:) < 0;
        h{i}.CData(switch_idx,:) = h{i}.CData(switch_idx,:)+repmat([0 0 1],sum(switch_idx),1);
    end
end

end


function [iScoreLocs,delta] = iScore(H,iBoxSize,deltaWindow,plotFlag)
%iScore computes insulation score from the Hi-C matrix
%   This function computes "insulation score" as described in Crane, Emily,
%   et al. "Condensin-driven remodelling of X chromosome topology during
%   dosage compensation." Nature 523.7559 (2015): 240. Insulation score can
%   be used to define TAD boundaries.
%   
%   H: raw Hi-C matrix, intra-chr (1D normalized, e.g. KR/VC/ICE, NOT O/E)
%   
%   iBoxSize: box size for I score (recommended 500kb, 5)
%   
%   deltaWindow: window size for delta score (recommended 100kb left/
%   right, 1)
%   
%   Scott Ronquist, 3/27/18

%% set default parameters (assume 100kb input data)
if nargin<4; plotFlag=1; end
if nargin<3; deltaWindow=1; end
if nargin<2; iBoxSize=5; end

%% i-score profile
% make sure symmetric
H = max(cat(3,H,H),[],3);

% calculate i-score
i_score = nan(length(H),1);
count = 1;
for i = (1+iBoxSize):(length(H)-(iBoxSize))
    box_locs_row = count:count+(iBoxSize-1);
    box_locs_col = i+1:(i+1)+(iBoxSize-1);
    temp = H(box_locs_row,box_locs_col);
    i_score(i) = nanmean(temp(:));
    count = count+1;
end

%% delta
i_score_norm = log2(i_score/nanmean(i_score));
delta = nan(length(i_score),1);
for i = ((1+iBoxSize)+1):((length(H)-(iBoxSize))-1)
    delta(i) = nanmean(i_score_norm(i-deltaWindow:i-1))-...
        nanmean(i_score_norm(i+1:i+deltaWindow));
end

% we want to find where delta is crossing zero, going down (i_score local minima)
iScoreLocs = find(diff(delta>0)<0);

% boundary strength
[delta_local_max,delta_local_max_loc] = findpeaks(delta);
[delta_local_min,delta_local_min_loc] = findpeaks(-delta);
delta_local_min = -delta_local_min;

boundary_strength = zeros(length(iScoreLocs),1);
for i = 1:length(iScoreLocs)-1
    temp = -(delta_local_max_loc-iScoreLocs(i));
    temp_max_loc = find(delta_local_max_loc==iScoreLocs(i)-min(temp(temp > 0)));
    
    temp = delta_local_min_loc-iScoreLocs(i);
    temp_min_loc = find(delta_local_min_loc==iScoreLocs(i)+min(temp(temp > 0)));
    
    boundary_strength(i) = delta_local_max(temp_max_loc)-...
        delta_local_min(temp_min_loc);
end

% filter
boundary_strength_thresh = .1;
iScoreLocs(boundary_strength<boundary_strength_thresh) = [];

%% plot (optional)
if 1==plotFlag
    figure, imagesc(log(H)), hold on
    for ii = 1:length(iScoreLocs)-1
        plot([iScoreLocs(ii)-.5 iScoreLocs(ii+1)-.5],[iScoreLocs(ii)-.5 iScoreLocs(ii)-.5], 'r-')
        plot([iScoreLocs(ii)-.5 iScoreLocs(ii+1)-.5],[iScoreLocs(ii+1)-.5 iScoreLocs(ii+1)-.5], 'r-')
        plot([iScoreLocs(ii)-.5 iScoreLocs(ii)-.5],[iScoreLocs(ii)-.5 iScoreLocs(ii+1)-.5], 'r-')
        plot([iScoreLocs(ii+1)-.5 iScoreLocs(ii+1)-.5],[iScoreLocs(ii)-.5 iScoreLocs(ii+1)-.5], 'r-')
    end
end

end


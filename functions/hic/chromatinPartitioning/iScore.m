function [iScoreLocs,delta] = iScore(H,iBoxSize,deltaWindow,plotFlag)
%iScore computes insulation score from the Hi-C matrix
%   This function computes "insulation score" as described in Crane, Emily,
%   et al. "Condensin-driven remodelling of X chromosome topology during
%   dosage compensation." Nature 523.7559 (2015): 240. Insulation score can
%   be used to define TAD boundaries.
%   
%   Input
%   H:              Raw Hi-C matrix, intra-chr (1D normalized, e.g. KR/VC/ICE, NOT O/E)
%   iBoxSize:       Box size for I score (recommended 500kb, 5)
%   deltaWindow:    Window size for delta score (recommended 100kb left/
%                   right, 1)
%   plotFlag:       Logical ([0 1]) to set whether to plot findings
%
%   Output
%   iScoreLocs:     Location of I-score boundaries
%   delta:          Vector of delta scores
%   
%   Scott Ronquist, 1/22/19

%% set default parameters (assume 100kb input data)
if nargin<4; plotFlag=1; end
if nargin<3; deltaWindow=1; end
if nargin<2; iBoxSize=5; end

%% i-score profile
% make sure symmetric
H = max(cat(3,H,H),[],3);

% calculate i-score
iScore = nan(length(H),1);
count = 1;
for i = (1+iBoxSize):(length(H)-(iBoxSize))
    boxLocsRow = count:count+(iBoxSize-1);
    boxLocsCol = i+1:(i+1)+(iBoxSize-1);
    temp = H(boxLocsRow,boxLocsCol);
    iScore(i) = nanmean(temp(:));
    count = count+1;
end

%% delta
iScoreNorm = log2(iScore/nanmean(iScore));
delta = nan(length(iScore),1);
for i = ((1+iBoxSize)+1):((length(H)-(iBoxSize))-1)
    delta(i) = nanmean(iScoreNorm(i-deltaWindow:i-1))-...
        nanmean(iScoreNorm(i+1:i+deltaWindow));
end

% we want to find where delta is crossing zero, going down (i_score local minima)
iScoreLocs = find(diff(delta>0)<0);

% remove start and end
iScoreLocs = iScoreLocs(2:end-1);

% boundary strength
[deltaLocalMax,deltaLocalMaxLoc] = findpeaks(delta);
[deltaLocalMin,deltaLocalMinLoc] = findpeaks(-delta);
deltaLocalMin = -deltaLocalMin;

boundaryStrength = zeros(length(iScoreLocs),1);
for i = 1:length(iScoreLocs)-1
    
    % find maxLoc left of iScoreLocs(i)
    temp = -(deltaLocalMaxLoc-iScoreLocs(i));
    tempMaxLoc = find(deltaLocalMaxLoc==iScoreLocs(i)-min(temp(temp >= 0)));
    
    % find minLoc right of iScoreLocs(i)
    temp = deltaLocalMinLoc-iScoreLocs(i);
    tempMinLoc = find(deltaLocalMinLoc==iScoreLocs(i)+min(temp(temp >= 0)));
    
    % determine boundary strength
    boundaryStrength(i) = deltaLocalMax(tempMaxLoc)-...
        deltaLocalMin(tempMinLoc);
end

% filter
boundaryStrengthThresh = .1;
iScoreLocs(boundaryStrength<boundaryStrengthThresh) = [];

%% plot (optional)
if 1==plotFlag
    figure, imagesc(log(H)), hold on
    for ii = 1:length(iScoreLocs)-1
        plot([iScoreLocs(ii)+.5 iScoreLocs(ii+1)+.5],[iScoreLocs(ii)+.5 iScoreLocs(ii)+.5], 'r-')
        plot([iScoreLocs(ii)+.5 iScoreLocs(ii+1)+.5],[iScoreLocs(ii+1)+.5 iScoreLocs(ii+1)+.5], 'r-')
        plot([iScoreLocs(ii)+.5 iScoreLocs(ii)+.5],[iScoreLocs(ii)+.5 iScoreLocs(ii+1)+.5], 'r-')
        plot([iScoreLocs(ii+1)+.5 iScoreLocs(ii+1)+.5],[iScoreLocs(ii)+.5 iScoreLocs(ii+1)+.5], 'r-')
    end
end

end


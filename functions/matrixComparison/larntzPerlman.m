function [H0,P,SMat] = larntzPerlman(R,n,alphaParam,plotFlag)
%larntzPerlman Larntz-Perlman procedure for testing covariance matrix equivalence
%
%   Input
%   R:          Sample covariance matrices, p x p x k
%   n:          Sample size
%   alphaParam: Alpha parameter to determine significance
%   plotFlag:   Logical for plotting max(Sij) vs chi-squared distribution
%   
%   Output
%   H0:         logical on whether the null hypothesis was rejected.
%               0 = same, 1 = different
%   P:          matrix with p-values, testing Sij in chi-squared
%               distribution (depreciated, used for testing purposes)
%   SMat:       S matrix
%
%   Example:
%   clear,rng(1)
%   n = 20;
%   p = 10;
%   a = rand(n,p);A = cov(a);
%   b = rand(n,p);B = cov(b);
%   [H0,P,SMat] = larntzPerlman(cat(3,A,B),n,.95,1);
%   
%   c = randn(n,p);C = cov(c);
%   [H0,P,SMat] = larntzPerlman(cat(3,A,B,C),n,.95,1);
%
%   Reference: Koziol, James A., et al. "A graphical technique for
%   displaying correlation matrices." The American Statistician 51.4
%   (1997): 301-304.
%   https://amstat.tandfonline.com/doi/pdf/10.1080/00031305.1997.10474402?needAccess=true
%
%   "if homogeneity is not rejected, data can be pooled"
%
%   Scott Ronquist, 1/22/19. scotronq@umich.edu

%% set default parameters
if nargin < 3;alphaParam = .95;fprintf('default alpha=%.2f...\n',alphaParam),pause(1);end
if nargin < 4;plotFlag = 0;fprintf('default plotFlag=%i...\n',plotFlag),pause(1);end

%% get upper triangle of matrices to vector
p = size(R,1);
k = size(R,3);

zLength = nchoosek(p,2);
z = zeros(zLength,k);

for iK = 1:k
    Rk = R(:,:,iK);
    mask = triu(true(size(Rk)),1);
    r = Rk(mask);
    z(:,iK) = .5*log((1+r)./(1-r));
    
    % Anderson-Darling test for normality % figure, histogram(z(:,iK))
    [adtestH,adtestP] = adtest(z(:,iK)); % 1 = normality rejected
end

zBar = mean(z,2);
SVec = (n-3)*sum((z-repmat(zBar,1,k)).^2,2);
T = max(SVec);
eAlpha = (1-alphaParam).^(2/(p*(p-1)));

%% test null hypothesis
 % https://www.mathworks.com/matlabcentral/answers/74472-how-can-i-create-the-chisquare-table
H0 = T > chi2inv(eAlpha,k-1);

% plot location in Chi-squared distribution
if plotFlag
    chiTemp = chi2pdf(.1:.1:T+5,k-1);
    figure, plot(.1:.1:T+5,chiTemp), hold on
    plot([T T], [0 nanmax(chiTemp)],'r-')
end

% P-Value matrix matrix
P = zeros(p, p);
P(triu(true(p),1)) = 1-chi2cdf(SVec,k-1);
P = P + P';

% reconstruct S as a matrix
SMat = zeros(p, p);
SMat(triu(true(p),1)) = SVec;
SMat = SMat + SMat';

end


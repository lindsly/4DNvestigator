function [H0,P,SMat] = larntzPerlman(R,n,alphaParam,plotFlag)
%larntzPerlman Larntz-Perlman procedure for testing correlation matrix equivalence
%
%   Input
%   R:          Sample correlation matrices, p x p x k
%   n:          Sample size
%   alphaParam: Alpha parameter to determine significance
%   plotFlag:   Logical for plotting max(Sij) vs chi-squared distribution
%   
%   Output
%   H0:         logical on whether the null hypothesis was rejected.
%               H_0: R^(1)=...R^(k)
%               0 = false, reject, data cannot be pooled
%               1 = true, accept, data can be pooled
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
%   Version 1.1 (4/26/19)
%   Written by: Scott Ronquist
%   Contact:    scotronq@umich.edu
%   Created:    1/22/19
%   
%   Revision History:
%   v1.0 (1/22/19)
%   * larntzPerlman.m created
%   v1.1 (4/26/19)
%   * update: code commenting

%% Set default parameters
% default alpha set to .95. No chi-squared plot by default
if nargin < 2
    n = size(R,1);
    fprintf('Default n=%i (size of input correlation matrices...\n',n)
    fprintf('Note: this default is specific for Hi-C matrices\n')
    pause(1);
end
if nargin < 3; alphaParam = .95; fprintf('Default alpha=%.2f...\n',alphaParam), pause(1); end
if nargin < 4; plotFlag = 0; fprintf('Default plotFlag=%i...\n',plotFlag), pause(1); end

%% Perform Larntz-Perlman procedure
% Get input data dimensions
p = size(R,1);
k = size(R,3);

% Fisher Z transform
zLength = nchoosek(p,2);
z = zeros(zLength,k);
for iK = 1:k
    Rk = R(:,:,iK);
    mask = triu(true(size(Rk)),1);
    r = Rk(mask);
    z(:,iK) = .5*log((1+r)./(1-r));
end

% Calculate mean Z score, S, and Test statistic, and Sidak correction
zBar = mean(z,2);
SVec = (n-3)*sum((z-repmat(zBar,1,k)).^2,2);
T = max(SVec);
eAlpha = (1-alphaParam).^(2/(p*(p-1)));

%% Test null hypothesis
% Chi-square inverse cumulative distribution function
fprintf('H_0 is rejected at level \x3B1 if T > X^2_((k-1),(\x3B5(\x3B1))\n')
if T > chi2inv(eAlpha,k-1)
    H0 = 0;
    fprintf('%.4f > %.4f, therefore H_0 = %i\n',T,chi2inv(eAlpha,k-1),H0)
else
    H0 = 1;
    fprintf('%.4f < %.4f, therefore H_0 = %i\n',T,chi2inv(eAlpha,k-1),H0)
end

% plot location in Chi-squared distribution
if plotFlag
    chiTemp = chi2pdf(.1:.1:T+5,k-1);
    figure, plot(.1:.1:T+5,chiTemp), hold on
    plot([T T], [0 nanmax(chiTemp)],'r-')
end

%% Format output
% P-Value matrix matrix
P = zeros(p, p);
P(triu(true(p),1)) = 1-chi2cdf(SVec,k-1);
P = P + P';

% reconstruct S as a matrix
SMat = zeros(p, p);
SMat(triu(true(p),1)) = SVec;
SMat = SMat + SMat';

end


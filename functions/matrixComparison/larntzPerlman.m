function [H0,M] = larntzPerlman(R,n,alphaParam)
%larntzPerlman Larntz-Perlman procedure for covariance matrix equivalence
%   R: sample covariance matrices, p x p x k
%   n: sample size
%   alphaParam: alpha parameter to determine significance
%   
%   output
%   H0: logical on whether the null hypothesis was rejected.
%   0 = same, 1 = different
%
%   Z: output z-transform matrix, p x p x k
%
%   Example:
%   clear,rng(1)
%   a = rand(10,5);A = cov(a);
%   b = rand(10,5);B = cov(b);
%   D = larntzPerlman(cat(3,A,B),10);
%
%   Reference: Koziol, James A., et al. "A graphical technique for
%   displaying correlation matrices." The American Statistician 51.4
%   (1997): 301-304.
%
%   "if homogeneity is not rejected, data can be pooled"
%
%   Scott Ronquist, 6/26/18

%% cat if nargin < 2
if nargin < 2;n = 100;fprintf('sample size not defined...\n'),pause(1);end
if nargin < 3;alphaParam = .95;fprintf('default alpha=%.2f...\n',alphaParam),pause(1);end

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
end

zBar = mean(z,2);
S = (n-3)*sum((z-repmat(zBar,1,k)).^2,2);
T = max(S);
eAlpha = (1-alphaParam).^(2/(p*(p-1)));

%% test null hypothesis
H0 = T > chi2inv(eAlpha,k-1); % https://www.mathworks.com/matlabcentral/answers/74472-how-can-i-create-the-chisquare-table

% alpha score matrix
M = zeros(p, p);
M(triu(true(p),1)) = chi2cdf(S,k-1);
M = M + M';

end


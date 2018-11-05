function H0 = larntzPerlman(A,n,alphaParam)
%larntzPerlman Larntz-Perlman procedure for covariance matrix equivalence
%   A: sample covariance matrices, p x p x k
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
%   if homogeneity is not rejected, data can be pooled
%
%   Scott Ronquist, 6/26/18

%% cat if nargin < 2
if nargin < 2;n = 100;fprintf('sample size not defined...\n'),pause(1);end
if nargin < 3;alphaParam = .95;fprintf('default alpha=%.2f...\n',alphaParam),pause(1);end

%% get upper triangle of matrices to vector
p = size(A,1);
k = size(A,3);

zLength = nchoosek(p,2);
z = zeros(zLength,k);

for iK = 1:k
    Ak = A(:,:,iK);
    mask = triu(true(size(Ak)),1);
    out = Ak(mask);
    z(:,iK) = .5*log((1+out)./(1-out));
end

zBar = mean(z,2);
S = (n-3)*sum((z-repmat(zBar,1,k)).^2,2);
T = max(S);
eAlpha = (1-alphaParam).^(2/(p*(p-1)));

H0 = T > chi2pdf(eAlpha,k-1); % this part still probably incorrect

end


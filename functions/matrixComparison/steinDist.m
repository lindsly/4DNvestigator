function [D] = steinDist(A,B)
%stein_dist This function calculates the Stein Distance between Matrices
%   A: sample covariance matrix, N x N
%   B: sample covariance matrix, N x N
%
%   Example:
%   rng(1)
%   a = rand(10,5);A = cov(a);
%   b = rand(10,5);B = cov(b);
%   D = steinDist(A,B);
%
%   Reference:
%   https://mail.google.com/mail/u/1/#search/stein/1642e6f30523c278
%
%   Scott Ronquist, 6/26/18

% D = trace((A-B)*pinv(A)*(A-B)*pinv(B));
D = trace(A*pinv(B))+trace(pinv(A)*B)-2*size(A,1);

end


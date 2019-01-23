function [D] = steinDist(A,B)
%steinDist This function calculates the Stein Distance between Matrices
%
%   Input
%   A:  Sample covariance matrix, N x N
%   B:  Sample covariance matrix, N x N
%
%   Output
%   D:  Stein Distance
%
%   Example:
%   rng(1)
%   a = rand(10,5);A = cov(a);
%   b = rand(10,5);B = cov(b);
%   D = steinDist(A,B);
%
%   Scott Ronquist, 1/22/19

%%
D = trace(A*pinv(B))+trace(pinv(A)*B)-2*size(A,1);

end


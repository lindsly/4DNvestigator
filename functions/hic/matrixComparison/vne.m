function [vneOut] = vne(A)
%vne Calculates the vonNeumann Entropy of A
%   
%   Inputs
%   A:          Input matrix	
%   
%   Outputs
%   vneOut:     Vonn Neumann entropy of A
%   
%   Version 1.0 (04/23/19)
%   Written by: Scott Ronquist
%   Contact: 	scotronq@umich.edu
%   Created: 	04/23/19
%   
%   Revision History:
%   v1.0 (04/23/19)
%   * vne.m created

%% set default parameters
% get laplacian
tempEig = eig(A);
tempEig = tempEig/(sum(tempEig));
vneOut = -sum(tempEig.*log(tempEig));

end

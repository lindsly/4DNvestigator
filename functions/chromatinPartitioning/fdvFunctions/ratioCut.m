% Function implementing the Ratio Cut algorithm which tends
% toward partitions of an equal or nearly equal size.
% RatioCut = Cut(A,B){1/(Vol(A)*Vol(B))}
% where RatioCut is minimized. Input is sortedI = I(p,p), the
% graph I sorted based on the sorted Fiedler vector

function [A,B] = ratioCut(sortedI)
numRows = size(sortedI,1);
cutVector = zeros(1,1);
breakpoint = floor(sqrt(numRows));
for j = breakpoint:numRows-breakpoint;
    cutAB = cut(sortedI,j,numRows);
    cutVector(j-breakpoint+1) = (cutAB/j) + (cutAB/(numRows-j));
end
[~, cutPoint] = min(cutVector);
cutPoint = cutPoint + breakpoint -1;
A = sortedI(1:cutPoint,1:cutPoint);
B = sortedI(cutPoint+1:numRows, cutPoint+1:numRows);


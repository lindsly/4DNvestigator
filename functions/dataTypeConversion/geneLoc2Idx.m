function [genePosNew] = geneLoc2Idx(genePos,binSize,trimLocs)
%geneLoc2Idx takes in gene position information and outputs trimmed Hi-C
%matrix indices
%
%   Input
%   genePos:    array with gene start and stop positions (Nx2 double)
%   binSize:    Hi-C bin size, in bp (scalar)
%   trimLocs:   Hi-C trimmed locations (Mx1 logical)
%
%   Output
%   genePosNew: array with gene start and stop positions, relative to
%               trimmed Hi-C matrix (Nx2 double)
%
%   Scott Ronquist, scotronq@umich.edu. 2/8/19

%% convert genePos bp -> Hi-C index
genePosNew = ceil(genePos/binSize);

%% shift positions with trimLocs
trimLocIdx = sort(find(trimLocs),'descend');

for iTrimLocs = 1:length(trimLocIdx)
    genePosNew(genePosNew >= trimLocIdx(iTrimLocs)) =...
        genePosNew(genePosNew >= trimLocIdx(iTrimLocs))-1;
end


end


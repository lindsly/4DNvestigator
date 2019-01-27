function [structOut] = parseHicHeaderInfo(cmdout)
%parseHicHeaderInfo parses the information output from read_hic_header.py
%
%   Input
%   cmdIn:      Characters output from read_hic_header.py
%
%   Output
%   structOut:  Structured information output
%
%   Scott Ronquist, scotronq@umich.edu. 1/22/19

%% Start
% split character string on new line character
C = strsplit(cmdout,newline);

% find which cells contain ":"
fieldIdx = [find(contains(C,':')),length(C)+1];

% create structure for output
structOut = struct;
for iField = 1:length(fieldIdx)-1
    colonLoc = strfind(C{fieldIdx(iField)},':');
    varName = C{fieldIdx(iField)}(1:colonLoc-1);
    varName = strrep(varName,'-','');
    varName = genvarname(varName);
    
    % format output
    switch varName
        case {'HiCVersion','MasterIndex'}
            structOut.(varName) =...
                str2double(C{fieldIdx(iField)+1:fieldIdx(iField+1)-1});
        case {'GenomeID'}
            fullPath = C{fieldIdx(iField)+1:fieldIdx(iField+1)-1};
            structOut.(varName) = fullPath;
            
            refGenomes = {'hg19','hg38','mm9','mm10'};
            for iRefGenomes = 1:length(refGenomes)
                if contains(fullPath,refGenomes{iRefGenomes})
                    structOut.refGenome = refGenomes{iRefGenomes};
                end
            end
            
        case {'Chromosomes'}
            chrsAvail = C(fieldIdx(iField)+1:fieldIdx(iField+1)-1);
            chr = {};
            chrLength = [];
            for iChr = 1:length(chrsAvail)
                chrsAvailSplit = strsplit(chrsAvail{iChr}, ' ');
                chrsAvailSplit(cellfun(@length,chrsAvailSplit)==0) = [];
                chr = [chr;chrsAvailSplit{1}];
                chrLength = [chrLength;str2num(chrsAvailSplit{2})];
            end
            chrSizes = table(chr, chrLength);
            
            structOut.(varName) = chrSizes;
            
        case {'BasePairdelimitedResolutions','FragmentdelimitedResolutions'}
            temp = str2double(C(fieldIdx(iField)+1:fieldIdx(iField+1)-1)');
            temp(isnan(temp)) = [];
            structOut.(varName) = temp;
            
        otherwise
            structOut.(varName) =...
                C(fieldIdx(iField)+1:fieldIdx(iField+1)-1)';
    end
    
end

end


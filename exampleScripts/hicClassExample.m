% This file will load all necessary information from the .hic file to the
% hic class
%
% Scott Ronquist, scotronq@umich.edu. 5/27/19

%% Set up
% Raw data (.hic and .results) files can be found here:
%   https://drive.google.com/open?id=1lSyU-7I0ME3X70Mt_-HjLMtPc-BMKHxm
clear
close all

% add path to Hi-C and RNA-seq data

%% Start
% Eventual input parameters for this function
%%%%%%%%%%%%%%%%%%%%%%%%%%%
juicer2matlabDir = fullfile('E:','MATLAB','4DNvestigator','functions','hic','juicer2matlab');
hicPath = 'Sample_64585_trim.hic';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize "hic" class variable
hicData = hic;

% Paths to juicer_tools scripts
juicerToolsPath = fullfile(juicer2matlabDir,'juicer_tools.jar');
juicerReadHicPath = fullfile(juicer2matlabDir,'straw-master','python','read_hic_header.py');

% Default parameters to extract
norm3d = {'observed','oe'};
norm1d = {'NONE','KR'};
binSize = 100E3;

% Temp files to output .hic data to
tempFnDump = [tempdir, 'tempDump.txt'];
tempFnNormKr = [tempdir, 'tempNormKr.txt'];
tempFnNormKrE = [tempdir, 'tempNormKrE.txt'];

%% Read .hic header
% Get hic header
[status,cmdout] = system(sprintf('python %s %s',juicerReadHicPath,hicPath));
if status
    error('Problem with reading hic header')
end

% Parse header info with 4DNvestigator function
[hicData.hicHeader] = parseHicHeaderInfo(cmdout);

%% Extract Hi-C header data
% Extract size close to 100kb
for iChr = 1:height(hicData.hicHeader.Chromosomes)
    tempChr = hicData.hicHeader.Chromosomes.chr{iChr};
    
    % Inter- vs intra-chr bp resolution default
    if strcmp(tempChr,'ALL')
        binSize = 1E6;
    else
        binSize = 100E3;
    end
    fprintf('Extracting data for chr:%s, resolution:%i...\n',...
        tempChr,binSize)
    
    %% Extract data usering juicer_tools
    % Extract raw data
    [statusDump,cmdoutDump] = system(sprintf('java -jar "%s" dump %s %s %s %s %s %s %s %s',...
        juicerToolsPath,'observed','NONE',hicPath,tempChr,tempChr,...
        'BP',num2str(binSize),tempFnDump));
    
    % Extract KR normalization vector
    [statusNormKr,cmdoutNormKr] = system(sprintf('java -jar "%s" dump %s %s %s %s %s %s %s %s',...
        juicerToolsPath,'norm','KR',hicPath,tempChr,...
        'BP',num2str(binSize),tempFnNormKr));
    
    % Extract E vector
    [statusNormKrE,cmdoutNormKrE] = system(sprintf('java -jar "%s" dump %s %s %s %s %s %s %s %s',...
        juicerToolsPath,'expected','KR',hicPath,tempChr,...
        'BP',num2str(binSize),tempFnNormKrE));
    
    %% Save to hic class
    % Read extracted data
    if statusDump~=0
        warning(cmdoutDump)
        hicData.rawData{end+1} = [];
    else
        tempRawData = readtable(tempFnDump,'Delimiter','\t'); delete(tempFnDump)
        if strcmp(tempChr,'ALL')
            hicData.rawData{end+1} = sparse(tempRawData{:,1}+1,tempRawData{:,2}+1,tempRawData{:,3});
        else
            hicData.rawData{end+1} = sparse((tempRawData{:,1}/binSize)+1,...
                (tempRawData{:,2}/binSize)+1,tempRawData{:,3});
        end
    end
    
    % Read extracted norm Kr
    if statusNormKr~=0
        warning(cmdoutNormKr)
        hicData.krVec{end+1} = [];
    else
        tempNormKr = readtable(tempFnNormKr,'Delimiter','\t'); delete(tempFnNormKr)
        if strcmp(tempChr,'ALL')
            hicData.krVec{end+1} = tempNormKr{2:end,1};
        else
            hicData.krVec{end+1} = tempNormKr{:,1};
        end
    end
    
    % Read extracted norm Kr E
    if statusNormKrE~=0
        warning(cmdoutNormKrE)
        hicData.krEVec{end+1} = [];
    else
        tempNormKrE = readtable(tempFnNormKrE,'Delimiter','\t'); delete(tempFnNormKrE)
        hicData.krEVec{end+1} = tempNormKrE{:,1};
    end
    
    % Add to data table
    hicData.rawDataTable = [hicData.rawDataTable; {tempChr,binSize}];
end

%% Extract data from Hi-C class
chrLoc1 = 21;
chrLoc2 = 21;
res = 100E3;
norm1d = 'KR';
norm3d = 'OE';

figure, imagesc(log(hicData.mat(chrLoc1,chrLoc2,res,norm1d,norm3d)))






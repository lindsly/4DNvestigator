% this file will load all necessary information from the .hic file to the
% hic class

clear
close all

%% Start
% eventual input parameters for this function
%%%%%%%%%%%%%%%%%%%%%%%%%%%
juicer2matlabDir = fullfile('E:','MATLAB','4DNvestigator','functions','hic','juicer2matlab');
hicPath = 'X:\projects\trisomy7_hcec\processed\hic\juicer\Sample_81952\inter_30.hic';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize "hic" class variable
hicData = hic;

% Paths to juicer_tools scripts
juicerToolsPath = fullfile(juicer2matlabDir,'juicer_tools.jar');
juicerReadHicPath = fullfile(juicer2matlabDir,'straw-master','python','read_hic_header.py');

% Default parameters to extract
norm3d = {'observed','oe'};
norm1d = {'NONE','KR'};
binSize = 100E3;

% temp files to load to
tempFnDump = [tempdir, 'tempDump.txt'];
tempFnNormKr = [tempdir, 'tempNormKr.txt'];
tempFnNormKrE = [tempdir, 'tempNormKrE.txt'];

%% Read .hic header
% Get hic header
[status,cmdout] = system(sprintf('python %s %s',juicerReadHicPath,hicPath));
if status
    error('Problem with reading hic header')
end

% Parse header info
[hicData.hicHeader] = parseHicHeaderInfo(cmdout);

%% Extract Hi-C header data
% Extract size close to 100kb
for iChr = 1:height(hicData.hicHeader.Chromosomes)
    tempChr = hicData.hicHeader.Chromosomes.chr{iChr};
    
    % inter vs intra resolution default
    if strcmp(tempChr,'ALL')
        binSize = 1E6;
    else
        binSize = 100E3;
    end
    fprintf('Extracting data for chr:%s, resolution:%i...\n',...
        tempChr,binSize)
    
    %% Extract data usering juicer_tools
    % Extract raw data
    [status,cmdout] = system(sprintf('java -jar "%s" dump %s %s %s %s %s %s %s %s',...
        juicerToolsPath,'observed','NONE',hicPath,tempChr,tempChr,...
        'BP',num2str(binSize),tempFnDump));
    if status~=0
        error(cmdout)
    end
    
    % Extract KR normalization vector
    [status,cmdout] = system(sprintf('java -jar "%s" dump %s %s %s %s %s %s %s %s',...
        juicerToolsPath,'norm','KR',hicPath,tempChr,...
        'BP',num2str(binSize),tempFnNormKr));
    if status~=0
        error(cmdout)
    end
    
    % Extract E vector
    [status,cmdout] = system(sprintf('java -jar "%s" dump %s %s %s %s %s %s %s %s',...
        juicerToolsPath,'expected','KR',hicPath,tempChr,...
        'BP',num2str(binSize),tempFnNormKrE));
    if status~=0
        error(cmdout)
    end
    
    %% Save to hic class
    % Read extracted data
    tempRawData = readtable(tempFnDump,'Delimiter','\t'); delete(tempFnDump)
    if strcmp(tempChr,'ALL')
        hicData.rawData{end+1} = sparse(tempRawData{:,1}+1,tempRawData{:,2}+1,tempRawData{:,3});
    else
        hicData.rawData{end+1} = sparse((tempRawData{:,1}/binSize)+1,...
            (tempRawData{:,2}/binSize)+1,tempRawData{:,3});
    end
    
    % Read extracted norm Kr
    tempNormKr = readtable(tempFnNormKr,'Delimiter','\t'); delete(tempFnNormKr)
    if strcmp(tempChr,'ALL')
        hicData.krVec{end+1} = tempNormKr{2:end,1};
    else
        hicData.krVec{end+1} = tempNormKr{:,1};
    end
    
    % Read extracted norm Kr E
    tempNormKrE = readtable(tempFnNormKrE,'Delimiter','\t'); delete(tempFnNormKrE)
    hicData.krEVec{end+1} = tempNormKrE{:,1};
    
    % Add to data table
    hicData.rawDataTable = [hicData.rawDataTable; {tempChr,binSize}];
end

















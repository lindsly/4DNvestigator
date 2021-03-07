
function [H] = fdnLoadGenericHic(Data_Loc,File_Name,iChr,binSize,iSample,H)
% fdnLoadGenericHic Loads and formats Hi-C data specified in dataInfo
% that is in the form: loc1 loc2 score
%
%   Input
%   Data_Loc:       string with folder where data is stored
%   File_Name:      string with file name
%   iChr:           current chromosome for this input file
%   binSize:        100 kb (1E5) or 1 Mb (1E6)
%   iSample:        current sample or time point
%   H:              Hi-C data structure where formatted data will be stored
%
%   Output
%   H:          structure containing all Hi-C data needed for 4DNvestigator
%
%   Written by: Stephen Lindsly
%   Contact:    lindsly@umich.edu
%   Created:    3/2/21
%   
%   Revision History:
%   v1.0 (3/7/21)
%   * fdnLoadGenericHic.m version 1 created


    %% default parameters
    if ~exist('numericType','var')||isempty(numericType);numericType='single';end
    if ~exist('H','var')||isempty(H);
        numChr = 24;
        H = fdnEmptyHMatrix(numChr);
    end

    % Get data from input file
    data = importdata([Data_Loc,File_Name]);
    intraFlag = 1;

    % Extract data from loc1, loc2, score formatted input
    H_temp = juicerDump2mat(data,intraFlag);
    
    % Trim data (removing rows/colums with mostly zeros) for better
    % processing later
    [H_temp, badLocs] = hicTrim(H_temp);

    % KR normalization and rescale back to original level (Note: may have
    % slightly different values than automatic KR normalization from Juicer)
    X = bnewt(H_temp);
    scale_factor1 = max(max(H_temp));
    tempKr = diag(X)*H_temp*diag(X);
    scale_factor2 = max(max(tempKr));
    scale_factor = scale_factor1/scale_factor2;
    tempKr(isnan(tempKr)) = 0;
    tempKr = tempKr*scale_factor;
    tempKr = max(cat(3,tempKr,tempKr'),[],3);

    % Observed/Expected normalization (Note: may have slightly different
    % values than automatic O/E normalization from Juicer)
    tempKrOe = ToepNorm(tempKr);
    tempKrOe(isnan(tempKrOe)) = 0;
    tempKrOe = max(cat(3,tempKrOe,tempKrOe'),[],3);
    tempKrE = tempKr./tempKrOe;

    % fix NaN locs derived from zero contacts in chr 1
    tempKrEVec = zeros(length(tempKrE),1);
    for iBin = 1:length(tempKrE)
        tempKrEVec(iBin) = nanmax(diag(tempKrE,iBin-1));
    end

    % interpolate NaN values
    nanx = isnan(tempKrEVec);
    t = 1:numel(tempKrEVec);
    tempKrEVec(nanx) = interp1(t(~nanx), tempKrEVec(~nanx), t(nanx));

    % fix NaNs at margin
    nanx = isnan(tempKrEVec);
    tempKrEVec(nanx) = tempKrEVec(find(nanx==0, 1, 'last' ));
    tempKrE = toeplitz(tempKrEVec);

    tempOE = tempKr./tempKrE(1:length(tempKr),1:length(tempKr));

    % change data type
    switch numericType
        case 'single'
            tempKr = single(tempKr);
            tempOE = single(tempOE);
    end

    % Assign to appropriate location 
    switch binSize
        case 1E5
            H.s100kb.krTrim{iChr}(:,:,iSample) = tempKr;
            H.s100kb.oeTrim{iChr}(:,:,iSample) = tempOE;
            H.s100kb.oeTrimBadLocs{iChr} = badLocs;
        case 1E6
            H.s1mb.krTrim{iChr}(:,:,iSample) = tempKr;
            H.s1mb.oeTrim{iChr}(:,:,iSample) = tempOE;
            H.s1mb.oeTrimBadLocs{iChr} = badLocs;
    end

end


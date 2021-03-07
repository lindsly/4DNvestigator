
function [H] = fdnLoadGenericHic(Data_Loc,File_Name,iChr,binSize,iSample,H)
% fdnLoadGenericHic Loads and formats Hi-C data specified in dataInfo
% that is in the form: loc1 loc2 contact
%
%   Input
%   dataInfo:       structure containing all experiment metadata
%   numericType:    string specifying numeric type to save as ('double'
%                   [default] or 'single')
%
%   Output
%   H:          structure containing all Hi-C data needed for 4DNvestigator
%
%   Written by: Stephen Lindsly
%   Contact:    lindsly@umich.edu
%   Created:    3/2/21
%   
%   Revision History:
%   v1.0 (3/2/21)
%   * fdnLoadGenericHic.m created


    %% default parameters
    if ~exist('numericType','var')||isempty(numericType);numericType='single';end
    if ~exist('H','var')||isempty(H);
        numChr = 24;
        H = fdnEmptyHMatrix(numChr);
    end

    data = importdata([Data_Loc,File_Name]);
    intraFlag = 1;

    H_temp = hicTrim(juicerDump2mat(data,intraFlag));
    [H_temp, badLocs] = hicTrim(H_temp);

    % extract KR
    X = bnewt(H_temp);
    scale_factor1 = max(max(H_temp));
    tempKr = diag(X)*H_temp*diag(X);
    scale_factor2 = max(max(tempKr));
    scale_factor = scale_factor1/scale_factor2;
    tempKr(isnan(tempKr)) = 0;
    tempKr = tempKr*scale_factor;
    tempKr = max(cat(3,tempKr,tempKr'),[],3);

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

    switch binSize
        case 1E5
            H.s100kb.krTrim{iChr}(1:length(tempKr),1:length(tempKr),iSample) = tempKr;
            H.s100kb.oeTrim{iChr}(1:length(tempKr),1:length(tempKr),iSample) = tempOE;
            H.s100kb.oeTrimBadLocs{iChr} = badLocs;
        case 1E6
            H.s1mb.krTrim{iChr}(1:length(tempKr),1:length(tempKr),iSample) = tempKr;
            H.s1mb.oeTrim{iChr}(1:length(tempKr),1:length(tempKr),iSample) = tempOE;
            H.s1mb.oeTrimBadLocs{iChr} = badLocs;
    end

end


function [h] = fdnSave(dataInfo,H,R)
%fdnSave saves important 4DNvestigator data structures
%   
%   Inputs
%   dataInfo:	
%   H:          Hi-C data structure
%   R:          RNA-seq data structure
%   
%   Outputs
%   h:          logical flag for succesful saving
%   
%   Version 1.0 (03/15/19)
%   Written by: Scott Ronquist
%   Contact:    scotronq@umich.edu
%   Created:    03/15/19
%   
%   Revision History:
%   v1.0 (03/15/19)
%   * fdnSave.m created

%% set default parameters
if ~exist('dataInfo','var')||isempty(dataInfo);dataInfo=[];end
if ~exist('H','var')||isempty(H);H=[];end
if ~exist('R','var')||isempty(R);R=[];end

%% format and save
% create new structure "saved Data"
savedData.dataInfo = dataInfo;
savedData.H = H;
savedData.R = R;

% save the data
try
    save(sprintf('%s%sdata%s%sData.mat',dataInfo.path.output,...
        dataInfo.delim,dataInfo.delim,...
        dataInfo.projName),'savedData','-v7.3')
    h = 1;
catch
    h = 0;
end

end

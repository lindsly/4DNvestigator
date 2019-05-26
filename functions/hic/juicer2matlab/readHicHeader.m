function [hicHeader] = readHicHeader(fn, idxRefGenome)
%readHicHeader reads the .hic header information
%
%   Input
%   fn:         full path to .hic file
%
%   Output
%   hicHeader: .hic header information
%
%   Scott Ronquist, scotronq@umich.edu. 1/22/19

%% format juicer tools input info
% find 4DNvestigator juicer tools path
juicerJarDir = mfilename('fullpath');

juicerJarDirLevels = strfind(juicerJarDir,filesep);
juicerJarDir = [juicerJarDir(1:juicerJarDirLevels(end)),'juicer_tools.jar'];

%% get and format .hic header info
% get "read_hic_header.py" file location
juicerReadHic = [juicerJarDir(1:juicerJarDirLevels(end)),...
    fullfile(
    sprintf('straw-master%spython%sread_hic_header.py',filesep,filesep)];

% run "read_hic_header.py"
[status,cmdout] = system(sprintf('python %s %s',juicerReadHic,fn));

% check alternate methods to obtain ref genome info
if status~=0
    if ~isempty(idxRefGenome)
        hicHeader.refGenome = idxRefGenome;
        hicHeader.Chromosomes = readtable(sprintf('%s.chrom.sizes',hicHeader.refGenome),...
            'fileType','text');
    else
        prompt = 'No header detected, please input reference genome name';
        opts = {'hg19','hg38'};
        [indx,tf] = listdlg('PromptString',prompt,'SelectionMode','single',...
            'ListString',opts);
        
        % get ref genome information
        hicHeader.refGenome = opts{indx};
        hicHeader.Chromosomes = readtable(sprintf('%s.chrom.sizes',opts{indx}),...
            'fileType','text');
    end
    
    hicHeader.Chromosomes.Properties.VariableNames = {'chr','chrLength'};
    hicHeader.Chromosomes.chr = cellfun(@(x) x(4:end), hicHeader.Chromosomes.chr, 'un', 0);
    
else
    % parse header info
    [hicHeader] = parseHicHeaderInfo(cmdout);
end

end


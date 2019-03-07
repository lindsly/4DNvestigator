function [hallmarkGsea] = readMSigDBHallmark()
%readMSigDBHallmark reads from MSig Database and formats output to
%structure
%
%   Input
%   (none)
%   
%   Output
%   hallmarkGsea: structure with MSigDB Hallmark genes structure
%
%   Scott Ronquist, scotronq@umich.edu. 12/20/18

%% MSigDB Hallmark Gene Set Collection
fid = fopen('h.all.v6.2.symbols.gmt');
tline = fgetl(fid);

hallmarkGsea = [];
while ischar(tline)
    tline = fgetl(fid);
    try
        tlineSplit = strsplit(tline,'\t');
    catch
        break
    end
    hallmarkGsea.(tlineSplit{1}) = sort(tlineSplit(3:end)');
end
fclose(fid);

end


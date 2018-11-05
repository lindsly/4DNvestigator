function [juicer_out] = juicer_tools_eig_mat(norm_1d,fn,chr,bp_frag,bin_size,fn_out)
%juicer_tools_dump_mat MATLAB version of juicertools "dump" to read .hic
%files
%   example: X:\projects\Wicha_BCSC\processed\hic\juicer\topdir_Sample_69655

if nargin<6;fn_out = 'juicer_temp.txt';end

%% num2str if necessary
if isnumeric(chr);chr=num2str(chr);end
if isnumeric(bin_size);bin_size=num2str(bin_size);end

%% call juicer tools
[a,b] = system(sprintf(['java -jar \\\\172.17.109.24\\internal_4DN\\tools\\juicer\\juicer_tools.jar ',...
    'eigenvector -p %s %s %s %s %s %s %s %s'],norm_1d,fn,chr,bp_frag,bin_size,fn_out));


% if ~isempty(fn_out)
    temp = readtable(fn_out,'ReadVariableNames',false);
    juicer_out = temp{:,1};
    if strcmp(fn_out,'juicer_temp.txt')
        delete juicer_temp.txt
    end
% else
%     line_locs = find(ismember(b, char([10 13])));
%     line_locs = [0,line_locs(1:end-1)];
%     juicer_out = zeros(length(line_locs),3);
%     for i = 1:length(line_locs)-1
%         fprintf('formatting juicer: %.4f\n',i/(length(line_locs)-1))
%         juicer_out(i,:) = str2num(b(line_locs(i)+1:line_locs(i+1)-1));
%     end
%     juicer_out(:,1:2) = round(juicer_out(:,1:2)/str2num(bin_size));
% end

end


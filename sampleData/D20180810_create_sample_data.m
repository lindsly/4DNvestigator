clear
close all

%% define parameters
param.norm_1d = 'KR';
param.norm_3d = 'oe';
param.bin_type = 'BP';
param.bin_size = 1E5;
param.chr = 14;
param.fn = 'X:\projects\Fibroblast_ts\processed\hic\juicer_analysis_S34030\inter_30.hic';

%% create sample data
A = juicer2mat(juicer_tools_dump_mat(param.norm_3d,param.norm_1d,...
    param.fn,param.chr,param.chr,param.bin_type,param.bin_size));

A(isnan(A)) = 0;
A = max(cat(3,A,A'),[],3);
A = hicTrim(A);

save('fibChr14KrOe100kbTrim.mat','A')

%% EXTRA
% A = A-diag(diag(A));
% A = hicTrim(A);
% ATrim = A;
% ATrim(ATrim<prctile(ATrim(:),95)) = 0;
% figure, plot(graph(ATrim),'EdgeAlpha',.1,'markerSize',10)
% [ABcomp,abGroup] = hicABcomp(A,[],[],[],[],1,'sign');
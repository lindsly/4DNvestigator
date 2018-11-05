% This script performs all default analysis methods for the Genome
% Structure Function Analysis Toolbox
% https://github.com/scotronq/srMatlabFunctions

% Scott Ronquist, 8/20/18

clear
close all
restoredefaultpath
addpath(genpath('.'))

%% Set data paths
param.projName = inputdlg('Input Project Name:','Project Name');
param.rnaseqPath = uigetdir(pwd,'Select RNA-seq data folder:');
param.hicPath = uigetdir(pwd,'Select Hi-C data folder:');
param.outputPath = uigetdir(pwd,'Select output directory');

addpath(genpath(param.rnaseqPath))
addpath(genpath(param.hicPath))

%% organize samples
% create sample info table, define samples/replicates
listing=dir(param.hicPath);[listing(:).dataType] = deal('hic');
sampleInfo=struct2table(listing(3:end));
listing=dir(param.rnaseqPath);[listing(:).dataType] = deal('rnaseq');
sampleInfo=[sampleInfo;struct2table(listing(3:end))];
sample = inputdlg(sampleInfo.name,'define samples',[1 50]);

% create replicate identifier
[uvals, ~, uidx] = unique(strcat(sample,sampleInfo.dataType));
replicate = ones(height(sampleInfo),1);    %mostly to copy the class and size
for K = 1 : length(uvals)
  mask = uidx == K;
  replicate(mask) = [1:sum(mask)]';
end

% create uniqname from info
uniqname = strcat(sampleInfo.dataType,{'_s'},sample,{'_r'},cellstr(num2str(replicate)));

sampleInfo = [sampleInfo,table(sample),table(replicate),table(uniqname)];

%% create data structure
fprintf('loading hic data...\n')
hicLocs = strcmp(sampleInfo.dataType,'hic');











data = cell2struct(cell(length(unique(sampleInfo.sample)),1),strcat({'s'},unique(sampleInfo.sample)),1);

fields = fieldnames(data);
for iSamp = 1:numel(fields)
    data.(fields{iSamp}) = struct('hic',[],'rnaseq',[]);
    
    iHicRes = [1E6,1E5];
    iHic1dNorm = {'KR','none'};
    iHic3dNorm = {'oe','observed'};
    for iHicNorm = 1:length(iHicRes)
        data.(fields{iSamp}).hic.b100kb.raw = hic2mat('observed','none',...
            [sampleInfo],'ALL','ALL','BP',1E5);
    end
    
end










param.dataResolution = [1E6,1E5];
param.dataNorm = {'none','observed';'KR','oe'};
uniqueSamples = unique(sampleInfo.sample);
for iChr = 1:22
    for iSamp = 1:length(uniqueSamples)
        for iNorm = 1:size(param.dataNorm,1)
            dataIdx = strcmp(sampleInfo.sample,'1') & strcmp(sampleInfo.dataType,'hic');
            
            data.hic(iChr,iSamp,iNorm) = data.(fields{iSamp}).hic.b100kb.raw = juicer2mat(...
            juicerToolsDump('observed','none',...
            [sample Info],...
            'ALL','ALL','BP',1E5));
        end
    end
end

data.hic(iChr,iSamp,iNorm)






param.bin_type = 'BP';
param.bin_size = [1E6,1E5];
param.norm_1d = 'KR';
param.norm_3d = 'oe';


%% Load data
hicLoc = find(strcmp(sampleInfo.sample,'hic'));
for iHic = 1:length(hicLoc)
    H.s1mb.KR = [];
    H.s1mb.KR_oe = [];
    for S = 1:length(samples.names)
        fprintf('extracting Hi-C 1MB genome-wide sample%d\n',S)
        H.s1mb.KR(:,:,S) = juicer2mat(...
            juicer_tools_dump_mat('observed','KR',...
            [selpath.hic,'/',samples.hic{S},'/inter_30.hic'],...
            'ALL','ALL','BP',1E6));
        
        temp = H.s1mb.KR(:,:,S);
        H.s1mb.KR_oe(:,:,S) = temp/mean(temp(temp>0));
        
        fprintf('chr:')
        for chr = 1:22%%NEED TO ADD CHRX/Y/M
            fprintf(' %d',chr)
            temp = H.s1mb.chr_info{chr,[2:3]};
            
            temp2 = juicer2mat(juicer_tools_dump_mat('observed','KR',...
                [selpath.hic,'/',samples.hic{S},'/inter_30.hic'],...
                chr,chr,'BP',1E6));
            
            H.s1mb.KR(temp(1):temp(1)+length(temp2)-1,...
                temp(1):temp(1)+length(temp2)-1,S) = temp2;
            
            temp2 = juicer2mat(juicer_tools_dump_mat('oe','KR',...
                [selpath.hic,'/',samples.hic{S},'/inter_30.hic'],...
                chr,chr,'BP',1E6));
            
            H.s1mb.KR_oe(temp(1):temp(1)+length(temp2)-1,...
                temp(1):temp(1)+length(temp2)-1,S) = temp2;
        end
        fprintf('\n')
    end
    H.s1mb.KR(isnan(H.s1mb.KR)) = 0;
    H.s1mb.KR_oe(isnan(H.s1mb.KR_oe)) = 0;
end
temp = juicer2mat(juicer_tools_dump_mat('observed',param.norm_1d,...
    [selpath.hic,'/',samples.hic{S},'/inter_30.hic'],Chr,Chr,param.bin_type,hic_param.bin_size));


%% Load RNA-seq data










%% Format Data
% this section formats the input files from Hi-C and RNA-seq pipelines and
% converts them to MATLAB files. All Hi-C data and information is contained
% within the "H" structure

H.s1mb.mat = [];
H.s100kb.mat = [];

for i = 1:length(samples.hic)
    % 1mb
    fprintf('loading Hi-C sample %s...\n',samples.hic{i})
    
    if ~isempty(strfind(samples.hic{i},'unknown'))
        continue
    end
    
    if exist([selpath.hic,'/',samples.hic{i},...
            '/extract_matrix-inter_30-observed-NONE-BP1MB/InterIntraChr.merged.txt'])~=2
        gunzip([selpath.hic,'/',samples.hic{i},...
            '/extract_matrix-inter_30-observed-NONE-BP1MB/InterIntraChr.merged.txt.gz']);
    end
    H.s1mb.chr_info = readtable([selpath.hic,'/',samples.hic{i},...
        '/extract_matrix-inter_30-observed-NONE-BP1MB/chrBinRangeInJuicer.txt']);
    H.s1mb.chr_info{:,2:3} = H.s1mb.chr_info{:,2:3}+1;
    H.s1mb.chr_info.Properties.VariableNames = {'chr','bin_start','bin_end'};
    H.s1mb.mat(:,:,i) = juicer2mat([selpath.hic,'/',samples.hic{i},...
        '/extract_matrix-inter_30-observed-NONE-BP1MB/InterIntraChr.merged.txt']);
    
    % 100kb
    if exist([selpath.hic,'/',samples.hic{i},...
            '/extract_matrix-inter_30-observed-NONE-BP100KB/InterIntraChr.merged.txt'])~=2
        gunzip([selpath.hic,'/',samples.hic{i},...
            '/extract_matrix-inter_30-observed-NONE-BP100KB/InterIntraChr.merged.txt.gz']);
    end
    H.s100kb.chr_info = readtable([selpath.hic,'/',samples.hic{i},...
        '/extract_matrix-inter_30-observed-NONE-BP100KB/chrBinRangeInJuicer.txt']);
    H.s100kb.chr_info{:,2:3} = H.s100kb.chr_info{:,2:3}+1;
    H.s100kb.chr_info.Properties.VariableNames = {'chr','bin_start','bin_end'};
end

%% RSEM matrix load
% this section loads the RNA-seq data output from RSEM. all RNA-seq data
% and information is contained within the "R" structure
fprintf('loading RNA-seq...\n')
R = rsem2mat(selpath.rnaseq,'hg19');
R_fields = fields(R);

for i = 1:length(R_fields)
    R.(R_fields{i})(~cellfun(@isnumeric,R.(R_fields{i}).chr),:) = [];
    R.(R_fields{i}).chr = cell2mat(R.(R_fields{i}).chr);
end

%100kb bins
chr_bin_lengths = diff(H.s100kb.chr_info{:,2:3},1,2)+1;
[R.s100kb.tpm,R.s100kb.gene] = rna2bin(R.TPM{:,7:end},R.TPM.gene_name,...
    [R.TPM.chr R.TPM.gene_start R.TPM.gene_end],1E5,chr_bin_lengths(1:22));

%1mb bins
chr_bin_lengths = diff(H.s1mb.chr_info{:,2:3},1,2)+1;
[R.s1mb.tpm,R.s1mb.gene] = rna2bin(R.TPM{:,7:end},R.TPM.gene_name,...
    [R.TPM.chr R.TPM.gene_start R.TPM.gene_end],1E6,chr_bin_lengths(1:22));

%% RNA-seq samples.hic aligned with Hi-C - WICHA SPECIFIC
R.TPM.Properties.VariableNames(7:end)'
samples.hic

if 1==1 % for Wicha samples
    samples.hic2RNAseq = {[7,8],[9,10],[4,5,6],[1,2,3]};
    samples.names = {'SPEN_control','SPEN_KD','ALDH_P','ALDH_N'};
else
    for i = 1:length(samples.hic)
        [samples.hic2RNAseq{i},tf] = listdlg('ListString',R.TPM.Properties.VariableNames(7:end)',...
            'PromptString',sprintf('Select RNA-seq for %s',samples.hic{i}),...
            'ListSize',[250,250]);
    end
end

fields = fieldnames(H);
for i = 1:numel(fields)
    for chr = 1:22
        temp = [];
        for S = 1:length(samples.hic)
            temp = [temp,...
                mean(R.(fields{i}).tpm{chr,1}(:,samples.hic2RNAseq{S}),2)];
        end
        R.(fields{i}).tpm_mean_sort{chr} = temp;
    end
end

%% Hi-C normalizations (juicer dump)
%%% PARAMETERS vvv
param.bin_type = 'BP';
param.bin_size = [1E6,1E5];
param.norm_1d = 'KR';
param.norm_3d = 'oe';
%%% PARAMETERS ^^^

fields = fieldnames(H);
for i = 1:numel(fields)
    H.(fields{i}).intra.juicer_oe = cell(22,1);
    for Chr = 1:22
        chr_size = diff(H.(fields{i}).chr_info{Chr,[2:3]})+1;
        H.(fields{i}).intra.juicer_oe{Chr} = zeros(chr_size,chr_size,length(samples.names));
        H.(fields{i}).intra.juicer_KR{Chr} = zeros(chr_size,chr_size,length(samples.names));
        
        for S = 1:length(samples.names)
            % KR
            % extract data from .hic to mat
            temp = juicer2mat(juicer_tools_dump_mat('observed',param.norm_1d,...
                [selpath.hic,'/',samples.hic{S},'/inter_30.hic'],Chr,Chr,param.bin_type,hic_param.bin_size));
            
            % concatenates and pads with NaN if different sizes
            H.(fields{i}).intra.juicer_KR{Chr} = padconcatenation_sr(...
                H.(fields{i}).intra.juicer_KR{Chr},temp,3);
            
            % KR and oe
            % extract data from .hic to mat
            temp = juicer2mat(juicer_tools_dump_mat(param.norm_3d,param.norm_1d,...
                [selpath.hic,'/',samples.hic{S},'/inter_30.hic'],Chr,Chr,param.bin_type,hic_param.bin_size));
            
            % concatenates and pads with NaN if different sizes
            H.(fields{i}).intra.juicer_oe{Chr} = padconcatenation_sr(...
                H.(fields{i}).intra.juicer_oe{Chr},temp,3);
        end
    end
end

%% Hi-C normalizations GENOME-WIDE (juicer dump)
if 1==0
    H.s1mb.KR = [];
    H.s1mb.KR_oe = [];
    for S = 1:length(samples.names)
        fprintf('extracting Hi-C 1MB genome-wide sample%d\n',S)
        H.s1mb.KR(:,:,S) = juicer2mat(...
            juicer_tools_dump_mat('observed','KR',...
            [selpath.hic,'/',samples.hic{S},'/inter_30.hic'],...
            'ALL','ALL','BP',1E6));
        
        temp = H.s1mb.KR(:,:,S);
        H.s1mb.KR_oe(:,:,S) = temp/mean(temp(temp>0));
        
        fprintf('chr:')
        for chr = 1:22%%NEED TO ADD CHRX/Y/M
            fprintf(' %d',chr)
            temp = H.s1mb.chr_info{chr,[2:3]};%diff(temp)+1
            
            temp2 = juicer2mat(juicer_tools_dump_mat('observed','KR',...
                [selpath.hic,'/',samples.hic{S},'/inter_30.hic'],...
                chr,chr,'BP',1E6));
            
            H.s1mb.KR(temp(1):temp(1)+length(temp2)-1,...
                temp(1):temp(1)+length(temp2)-1,S) = temp2;
            
            temp2 = juicer2mat(juicer_tools_dump_mat('oe','KR',...
                [selpath.hic,'/',samples.hic{S},'/inter_30.hic'],...
                chr,chr,'BP',1E6));
            
            H.s1mb.KR_oe(temp(1):temp(1)+length(temp2)-1,...
                temp(1):temp(1)+length(temp2)-1,S) = temp2;
        end
        fprintf('\n')
    end
    H.s1mb.KR(isnan(H.s1mb.KR)) = 0;
    H.s1mb.KR_oe(isnan(H.s1mb.KR_oe)) = 0;
end

%% Remove centromeres and low signal nodes
%%% PARAMETERS vvv
param.num_diag = 1;
param.num_sparse = 0;
%%% PARAMETERS ^^^

for i = 1:numel(fields)
    fprintf('normalizing Hi-C %s\n',fields{i})
    for Chr = 1:22
        % find bad locs
        [H_temp,bad_locs] = hicTrim(H.(fields{i}).intra.juicer_oe{Chr},...
            param.num_diag,param.num_sparse);
        
        % Trim all necessary data structures
        H.(fields{i}).intra.juicer_oe_trim{Chr} =...
            H.(fields{i}).intra.juicer_oe{Chr}(~bad_locs,~bad_locs,:);
        
        H.(fields{i}).intra.juicer_KR_trim{Chr} =...
            H.(fields{i}).intra.juicer_KR{Chr}(~bad_locs,~bad_locs,:);
        
        R.(fields{i}).tpm_juicer_trim{chr,1} =...
            R.(fields{i}).tpm{chr,1}(~bad_locs,:);
        
        R.(fields{i}).gene_juicer_trim{chr,1} =...
            R.(fields{i}).gene{chr,1}(~bad_locs,:);
        
        H.(fields{i}).intra.juicer_oe_zeros_locs{Chr} = bad_locs;
    end
end

%% Eigenspectrum VNGE and plots
%%% PARAMETERS vvv
param.num_eigs = 10;
param.H_norm_method = 'corr';
%%% PARAMETERS ^^^

eig_spec = zeros(length(samples.names),22,10);
VNE = zeros(length(samples.names),22);

fields = fieldnames(H);
for i = 1:numel(fields)
    for Sample = 1:length(samples.names)
        for Chr = 1:22
            fprintf('Calculating eigspec and VNE. Sample:%s (%d/%d), chr:%d...\n',...
                samples.names{Sample},Sample,length(samples.names),Chr)
            
            % calculate eig spec
            [VN_entropy, D_all] = hicVnEntropy(H.(fields{i}).intra.juicer_oe_trim{Chr}(:,:,Sample),...
                param.num_eigs, 1, param.H_norm_method);
            
            VNE(Sample,Chr) = VN_entropy;
            eig_spec(Sample,Chr,:) = D_all;
        end
    end
end

% figure eigspec
for Sample = 1:length(samples.names)   
    [h] = plot_hic_eigspec(squeeze(eig_spec(Sample,:,:)),samples.names{Sample});
end

% figure VNE compare
for Sample1 = 1:length(samples.names)
    for Sample2 = Sample1+1:length(samples.names)
        [h] = plot_hic_vne_diff([VNE(Sample1,:);VNE(Sample2,:)],...
            {samples.names{Sample1};samples.names{Sample2}});
    end
end

%% A/B compartments
fields = fieldnames(H);
for i = 1:numel(fields)
    for chr = 1:22
        for S = 1:length(samples.hic)
            fprintf('extracting AB compartments %s chr%d sample %d\n',fields{i},chr,S)
            [H.(fields{i}).intra.ab_comp_juicer{chr,S}] = hicABcomp(double(...
                H.(fields{i}).intra.juicer_oe_trim{chr}(:,:,S)),...
                'PC1',R.(fields{i}).tpm_juicer_trim_mean_sort{chr}(:,S));
            
            % normalize AB from -1 to 1
            if 1==1
                H.(fields{i}).intra.ab_comp_juicer{chr,S}=...
                    (H.(fields{i}).intra.ab_comp_juicer{chr,S}-...
                    min(H.(fields{i}).intra.ab_comp_juicer{chr,S}))/...
                    (max(H.(fields{i}).intra.ab_comp_juicer{chr,S})-...
                    min(H.(fields{i}).intra.ab_comp_juicer{chr,S}))*2-1;
            end
        end
        if 1==0
            figure,
            for S = 1:length(samples.hic)
                subplot(length(samples.hic),1,1),
                bar(H.(fields{i}).intra.ab_comp_juicer{chr,S})
                title(samples.names{S})
            end
            suptitle(sprintf('AB compartments,  chr: %d, resolution: %s',...
                chr,fields{i}))
        end
    end
end

%% whole genome AB figures (good figs in here)
if 1==0
    figure, plot(cell2mat(H.s100kb.intra.ab_comp_juicer)),legend(strrep(samples.names,'_',' '))
    figure, plot(cell2mat(H.s100kb.intra.ab_comp_juicer(9,:)))
    figure, imagesc(double(H.s1mb.intra.juicer_oe_trim{12}(:,:,3)))
    figure, plot(sum(double(H.s1mb.intra.juicer_oe_trim{12}(:,:,3))))
    
    % save subset of data for stephen
    juicer_oe=H.s100kb.intra.juicer_oe;
    juicer_oe_trim=H.s100kb.intra.juicer_oe_trim;
    scott_ab_pc1_scaled = H.s100kb.intra.ab_comp_juicer;
    save('wicha_bcsc_data_subset.mat','juicer_oe_trim','juicer_oe','scott_ab_pc1_scaled','R')
    
    % AB comparison between all samples
    fields = fieldnames(H);
    for i = 2%1:numel(fields)
        for chr = 1:22
            figure,
            for S = 1:length(samples.hic)
                subplot(length(samples.hic),1,S),
                bar(H.(fields{i}).intra.ab_comp_juicer{chr,S})
                title(strrep(samples.names{S},'_',' '))
            end
            linkaxes(get(gcf,'children'))
            suptitle(sprintf('AB compartments,  chr: %d, resolution: %s',...
                chr,fields{i}))
        end
    end
    
    % 1mb overview
    figure,
    for S = 1:length(samples.name)
        subplot(2,2,S),
        plot_hic(log(H.s1mb.mat(:,:,1)),'erez',0)
        %         axis square
        title(strrep(samples.names{S},'_',' '))
    end
    linkaxes(get(gcf,'children'))
    suptitle(sprintf('AB compartments,  chr: %d, resolution: %s',...
        chr,fields{i}))
    
    %% analysis example 1
    chr_select=14;
    sample_select=3;
    
    for sample_select=1:length(samples.names)
        figure('position',[100 100 700 1000]),
        subplot(6,1,1)
        temp = log2(R.s100kb.tpm_juicer_trim_mean_sort{chr_select}(:,sample_select)+.5)+1;
        bar(temp)
        xlim([0 length(temp)])
        ax = gca;ax.LineWidth = 2;ax.FontSize=8;
        ylabel('log2(TPM)','fontsize',12)
        title('RNA-seq','fontsize',15,'fontweight','bold')
        
        subplot(6,1,2)
        temp = H.s100kb.intra.ab_comp_juicer{chr_select,sample_select};
        bar(find(temp<0),temp(temp<0),'r'),hold on
        bar(find(temp>=0),temp(temp>=0),'g')
        xlim([0 length(temp)])
        ax = gca;ax.LineWidth = 2;ax.FontSize=8;
        ylabel('normalized [-1,1]','fontsize',12)
        title('Chromatin Patterning','fontsize',15,'fontweight','bold')
        
        subplot(6,1,3:6)
        temp = corr(corr(double(H.s100kb.intra.juicer_oe_trim{chr_select}(:,:,sample_select))));
        temp = temp-diag(diag(temp));
        imagesc(temp)
        ax = gca;ax.LineWidth = 2;ax.FontSize=8;
        title('Hi-C','fontsize',15,'fontweight','bold')
        
        ax = suptitle(sprintf('%s, chr%d, 100kb',...
            strrep(samples.names{sample_select},'_',' '),chr_select));
        ax.FontWeight='bold';
        saveas(gcf,[selpath.out,'\',sprintf('chr%d_sample%d.png',chr_select,sample_select)])
    end
    
    %% analysis example 2
    % all data together
    all_H = cell2mat(H.s100kb.intra.ab_comp_juicer);
    all_R = cell2mat(R.s100kb.tpm_juicer_trim_mean_sort');
    all_R = log(all_R+.5)+1;
    all_R = all_R*(mean(std(all_H))/mean(std(all_R)));
    all_gene = vertcat(R.s100kb.gene_juicer_trim{:});
    
    %spen control vs ALDH-
    s1=1;s2=4;
    dist_all_1 = zeros(size(all_R,1),1);
    for i=1:size(all_R,1)
        dist_all_1(i) = pdist([all_R(i,s1),all_H(i,s1);all_R(i,s2),all_H(i,s2)]);
    end
    sum(dist_all_1)
    
    %spen KD vs ALDH-
    s1=2;s2=4;
    dist_all_2 = zeros(size(all_R,1),1);
    for i=1:size(all_R,1)
        dist_all_2(i) = pdist([all_R(i,s1),all_H(i,s1);all_R(i,s2),all_H(i,s2)]);
    end
    sum(dist_all_2)
    
    [~,idx] = sort(dist_all_2-dist_all_1,'ascend');
    
    figure('position',[100 100 700 700]), hold on
    plot(all_R(idx(1:10),[1,2,4]),all_H(idx(1:10),[1,2,4]),'*')
    plot([all_R(idx(1:10),1),all_R(idx(1:10),2)]',...
        [all_H(idx(1:10),1),all_H(idx(1:10),2)]','k-')
    plot([all_R(idx(1:10),2),all_R(idx(1:10),4)]',...
        [all_H(idx(1:10),2),all_H(idx(1:10),4)]','k-')
    
    for i = 1:10
        text(mean([all_R(idx(i),1),all_R(idx(i),2)]),...
            mean([all_H(idx(i),1),all_H(idx(i),2)]),...
            all_gene{idx(i)})%,'backgroundcolor','w')
    end
    
    legend(strrep(samples.names([1,2,4]),'_',' '),'FontSize',14)
    axis equal
    
    %     ax = gca;ax.LineWidth = 2;ax.FontSize=8;
    ylabel('Chromatin Patterning, Normalized','fontweight','bold','fontsize',12)
    xlabel('RNA-seq, Normalized','fontweight','bold','fontsize',12)
    title('Structure Function change','fontsize',15,'fontweight','bold')
end

%% TAD boundaries
fields = fieldnames(H);
for i = 2%1:numel(fields)
    for chr = 1:22
        fprintf('extracting TAD boundaries, i-score %s chr%d\n',fields{i},chr)
        [H.(fields{i}).intra.iScore{chr},delta] = iScore(H.(fields{i}).intra.raw{chr});
        
        H.(fields{i}).intra.tad_laplace{chr} = TAD_Laplace_Sijia(H.(fields{i}).intra.toep_trim{chr});
    end
end

%% Centrality (UNDER CONSTRUCTION, works with test data)
test_data.size = 10;
test_data.TPs = 4;

test_data.hic = zeros(test_data.size,test_data.size,test_data.TPs);
test_data.rnaseq = zeros(test_data.size,test_data.TPs);

for i = 1:test_data.TPs
    temp = rand(test_data.size);
    test_data.hic(:,:,i) = temp'*temp;
    test_data.rnaseq(:,i) = rand(test_data.size,1);
end

[features] = sfAnalysis(test_data.hic,test_data.rnaseq);

%% SF analysis phaseplane UNDER CONSTRUCTION

fields = fieldnames(H);
for chr = 1:22
    fprintf('extracting SF phase plane %s chr%d\n',fields{i},chr)
    sf_phaseplane(R.s100kb.tpm_juicer_trim_mean_sort{chr},...
        cell2mat(H.s100kb.intra.ab_comp_juicer(chr,:)),...
        R.s100kb.gene_juicer_trim{chr},strrep(samples.names,'_',' '),10,'log2')
end

sf_phaseplane(R,H,bin_names,sample_names)

%% RNA-seq PCA
% log2 TPM: https://www.nature.com/articles/sdata2016109#f2
[coeff,score,latent,tsquared,explained] = pca(log2(R.TPM{:,7:end}'+.5)+1);
figure, hold on
for S = 1:length(samples.hic)
    scatter3(score(samples.hic2RNAseq{S},1),score(samples.hic2RNAseq{S},2),...
        score(samples.hic2RNAseq{S},3))
end
axis equal
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
zlabel('3rd Principal Component')
title('RNA-seq PCA, log2(TPM)')
legend(samples.names)

%% RNA-seq MATLAB DE (UNDER CONSTRUCTION, works with test data)
% https://www.mathworks.com/help/bioinfo/examples/identifying-differentially-expressed-genes-from-rna-seq-data.html

% like DEseq, negative binomial distribution, with variance and mean linked by local regression
% https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106

load pasilla_count_noMM.mat
counts = geneCountTable{:,3:end};
[pvalue,padj] = matlab_negbin_DE(counts(:,1:2),counts(:,3:4));


samples.names = strrep(samples.names,'+','P');
samples.names = strrep(samples.names,'-','N');
for i = 1:length(samples.names)
    for ii = i+1:length(samples.names)
        fprintf('running MATLAB DESeq for %s vs %s...\n',samples.names{i},samples.names{ii})
        [pvalue,padj] = matlab_negbin_DE(R.expected_count(:,6+samples.hic2RNAseq{i}),...
            R.expected_count(:,6+samples.hic2RNAseq{ii}));
        
        meanTreated = mean(R.TPM{:,6+samples.hic2RNAseq{i}},2);
        meanUntreated = mean(R.TPM{:,6+samples.hic2RNAseq{ii}},2);
        foldChange = meanTreated ./ meanUntreated;
        log2FC = log2(foldChange);
        gene_name = R.TPM.gene_name;
        chr = R.TPM.chr;
        gene_start = R.TPM.gene_start;
        gene_end = R.TPM.gene_end;
        
        R.matlabDE{i,ii}=table(gene_name,chr,gene_start,gene_end,...
            meanTreated,meanUntreated,foldChange,log2FC,pvalue,padj);
        
        R.matlabDE{i,ii}.Properties.VariableNames{5}=...
            samples.names{i}(find(~isspace(samples.names{i})));
        R.matlabDE{i,ii}.Properties.VariableNames{6}=...
            samples.names{ii}(find(~isspace(samples.names{ii})));
        
        R.matlabDE{i,ii} = sortrows(R.matlabDE{i,ii},'padj','ascend');
        
    end
end

%% GSEA
% Prior to conducting gene set enrichment analysis, conduct your
% differential expression analysis using any of the tools developed by the
% bioinformatics community (e.g., cuffdiff, edgeR, DESeq, etc).

% Based on your differential expression analysis, rank your features and
% capture your ranking in an RNK-formatted file. The ranking metric can be
% whatever measure of differential expression you choose from the output of
% your selected DE tool. For example, cuffdiff provides the (base 2) log of
% the fold change.

% Run GSEAPreranked, but make sure to select "classic" for your enrichment
% score (thus, not weighting each gene's contribution to the enrichment
% score by the value of its ranking metric).

for i = 1:length(samples.names)
    for ii = i+1:length(samples.names)
        
        temp = sortrows(R.matlabDE{i,ii},'log2FC','descend');
        temp(isnan(temp.log2FC),:) = [];
        temp(isinf(abs(temp.log2FC)),:) = [];
        temp(temp.padj>.05,:) = [];
        
        writetable(temp(:,[1,8]),[selpath.out,'\',sprintf('%s_over_%s.rnk',...
            temp.Properties.VariableNames{5},...
            temp.Properties.VariableNames{6})],...
            'WriteRowNames',false,'Delimiter','\t','FileType','text')
        
    end
end

%% HiCSpector
S_d = zeros(size(H.s100kb.mat,3),size(H.s100kb.mat,3),22);
Q = zeros(size(H.s100kb.mat,3),size(H.s100kb.mat,3),22);

for i = 1:size(H.s100kb.mat,3)
    for ii = i+1:size(H.s100kb.mat,3)
        for chr = 1:22
            [S_d(i,ii,chr),Q(i,ii,chr)] = hic_spector(double(H.s100kb.intra.juicer_oe{chr}(:,:,i)),...
                double(H.s100kb.intra.juicer_oe{chr}(:,:,i)));
        end
    end
end

%% SAVE output
if 1==0
    save([selpath.out,'\sf_pipeline_data.mat'],'H','R','samples','selpath','-v7.3')
end

%% Wicha Data (goog figs in here)
if 1==0
    % 69655	SUM159	ALDH+ (stem cell like)
    % 69656	SUM159	ALDH- (not stem cell like)
    % 103172	SUM159	SUM159-SPEN_ctr
    % 103173	SUM159	SUM159-SPEN-dox
    %%%%%%%%%%%%%%%%%
    samples.hic(:,2) = {'SPEN ctr';'SPEN KD';'ALDH+';'ALDH-'};
    parameters.num_eigs = 10;
    parameters.H_norm_method = 'corr';
    %%%%%%%%%%%%%%%%%
    
    H.(fields{i}).vn_entropy = zeros(height(H.(fields{i}).chr_info), length(samples.hic));
    H.(fields{i}).eig_spec = zeros(height(H.(fields{i}).chr_info), parameters.num_eigs,length(samples.hic));
    
    fields = fieldnames(H);
    for i = 1:numel(fields)
        fprintf('extracting eig-spectrum and VNGE Hi-C %s\n',fields{i})
        
        % calculate VNGE and eigenspectrum
        for chr = 1:22
            fprintf('extracting eig-spectrum and VNGE Hi-C %s chr%d\n',fields{i},chr)
            [H.(fields{i}).vn_entropy(chr,:), D_all] =...
                hicVnEntropy(double(H.(fields{i}).intra.juicer_oe_trim{chr}),...
                parameters.num_eigs, 1, parameters.H_norm_method);
            H.(fields{i}).eig_spec(chr,:,:) = reshape(D_all,1,parameters.num_eigs,length(samples.hic));
        end
        
        % figure
        figure
        title(sprintf('Hi-C %s VNGE per chr (juicer O/E, corr, 10 eigs)',fields{i}))
        bar(H.(fields{i}).vn_entropy)
        legend(samples.hic(:,2))
        
        figure
        title(sprintf('Hi-C %s VNGE total (juicer O/E, corr, 10 eigs)',fields{i}))
        bar(sum(H.(fields{i}).vn_entropy))
        xticklabels(samples.hic(:,2))
        
    end
    
    % figure
    chr = 2;
    bin_scale = 1;
    figure
    for S = 1:4
        subplot(2,2,S), imagesc(H.(fields{bin_scale}).intra.juicer_oe_trim{chr}(:,:,S))
        title(samples.hic(S,2))
    end
    linkaxes(get(gcf,'children'))
    suptitle(sprintf('chr: %d, scale: %s',chr,fields{bin_scale}))
    
    % whole genome VNGE, 1mb IN PROGRESS
    temp = min(sum(H.s1mb.mat),[],3)==0;
    temp2 = H.s1mb.mat;
    temp2(temp,:,:) = [];
    temp2(:,temp,:) = [];%figure, imagesc(log(temp2(:,:,1)))
    
    bin_info = [];%zeros(size(H.s1mb.mat,1),3);
    for i = 1:height(H.s1mb.chr_info)
        temp = [1:diff(H.s1mb.chr_info{i,2:3})+1]'*1E6;
        bin_info = [bin_info;[ones(length(temp),1)*i,temp-(1E6-1),temp]];
    end
    
    temp = norm_hic_bins(H.s1mb.mat,bin_info);%figure, imagesc(log(temp3(:,:,1)))
    
    temp2 = min(sum(temp),[],3)==0;
    temp3 = temp;
    temp3(temp2,:,:) = [];
    temp3(:,temp2,:) = [];
    
    [a,~] = hicVnEntropy(temp3(1:2733,1:2733,:),10,1,'corr');
    
    figure
    title(sprintf('Hi-C 1MB VNGE total (juicer, Toep O/E, corr, 20 eigs)'))
    bar(a)
    xticklabels(samples.hic(:,2))
end

%% remove centromeres and low signal nodes (OLD)
if 1==0
for i = 1:numel(fields)
    fprintf('normalizing Hi-C %s\n',fields{i})
    for chr = 1:22
        % trim Hi-C (centromere remove)
        
        temp = [];
        for S = 1:length(samples.names)
            num_obs = sum(logical(H.(fields{i}).intra.juicer_oe{chr}(:,:,S)-...
                diag(diag(H.(fields{i}).intra.juicer_oe{chr}(:,:,S)))-...
                diag(diag(H.(fields{i}).intra.juicer_oe{chr}(:,:,S),1),1)-...
                diag(diag(H.(fields{i}).intra.juicer_oe{chr}(:,:,S),-1),-1)));
            num_obs_low=num_obs<median(num_obs)-2*std(num_obs);
            num_obs_low_rel=num_obs<length(num_obs)*.05;
            
            num_counts = sum(H.(fields{i}).intra.juicer_oe{chr}(:,:,S)-...
                diag(diag(H.(fields{i}).intra.juicer_oe{chr}(:,:,S)))-...
                diag(diag(H.(fields{i}).intra.juicer_oe{chr}(:,:,S),1),1)-...
                diag(diag(H.(fields{i}).intra.juicer_oe{chr}(:,:,S),-1),-1));
            num_counts_zero = num_counts==0;
            num_counts_low = num_counts<mean(num_counts)-2*std(num_counts);
            num_counts_high = num_counts<mean(num_counts)-2*std(num_counts);
            
            temp = [temp;num_obs_low;num_obs_low_rel;...
                num_counts_zero;num_counts_low;num_counts_high];
        end
        
        temp2 = logical(sum(temp,1));
        H.(fields{i}).intra.juicer_oe_zeros_locs{chr} = temp2;
        H.(fields{i}).intra.juicer_oe_trim{chr} = H.(fields{i}).intra.juicer_oe{chr};
        H.(fields{i}).intra.juicer_oe_trim{chr}(temp2,:,:) = [];
        H.(fields{i}).intra.juicer_oe_trim{chr}(:,temp2,:) = [];
        
        H.(fields{i}).intra.juicer_KR_trim{chr} = H.(fields{i}).intra.juicer_KR{chr};
        H.(fields{i}).intra.juicer_KR_trim{chr}(temp2,:,:) = [];
        H.(fields{i}).intra.juicer_KR_trim{chr}(:,temp2,:) = [];
        
        % trim RNAseq (juicer)
        R.(fields{i}).tpm_juicer_trim{chr,1} = R.(fields{i}).tpm{chr,1};
        R.(fields{i}).tpm_juicer_trim{chr,1}(temp2,:) = [];
        
        R.(fields{i}).gene_juicer_trim{chr} = R.(fields{i}).gene{chr};
        R.(fields{i}).gene_juicer_trim{chr}(temp2) = [];
        
        % Wicha RNAseq mean sorted w/ samples.hic
        temp = [];
        for S = 1:length(samples.names)
            temp = [temp,...
                mean(R.(fields{i}).tpm_juicer_trim{chr,1}(:,samples.hic2RNAseq{S}),2)];
        end
        
        R.(fields{i}).tpm_juicer_trim_mean_sort{chr} = temp;
    end
end
end

%% Hi-C normalizations (juicer dump) - OLD
if 1==0
fields = fieldnames(H);
for i = 1:numel(fields)
    for chr = 1:22
        chr_size = diff(H.(fields{i}).chr_info{chr,[2:3]})+1;
        H.(fields{i}).intra.juicer_oe{chr} = zeros(chr_size,chr_size,length(samples.names));
        H.(fields{i}).intra.juicer_KR{chr} = zeros(chr_size,chr_size,length(samples.names));
        for S = 1:length(samples.names)
            fprintf('normalizing Hi-C %s chr%d sample%d\n',fields{i},chr,S)
            
            chr_start = H.(fields{i}).chr_info{chr,2};
            chr_end = H.(fields{i}).chr_info{chr,3};
            
            if strcmp(fields{i},'s1mb')
                temp = juicer2mat(...
                    juicer_tools_dump_mat('oe','KR',...
                    [selpath.hic,'/',samples.hic{S},'/inter_30.hic'],...
                    chr,chr,'BP',1E6));
                temp(isnan(temp)) = 0;
                H.(fields{i}).intra.juicer_oe{chr}(1:size(temp,1),1:size(temp,2),S) = temp;
                
                temp = juicer2mat(...
                    juicer_tools_dump_mat('observed','KR',...
                    [selpath.hic,'/',samples.hic{S},'/inter_30.hic'],...
                    chr,chr,'BP',1E6));
                temp(isnan(temp)) = 0;
                H.(fields{i}).intra.juicer_KR{chr}(1:size(temp,1),1:size(temp,2),S) = temp;
            else
                temp = juicer2mat(...
                    juicer_tools_dump_mat('oe','KR',...
                    [selpath.hic,'/',samples.hic{S},'/inter_30.hic'],...
                    chr,chr,'BP',1E5));
                temp(isnan(temp)) = 0;
                H.(fields{i}).intra.juicer_oe{chr}(1:size(temp,1),1:size(temp,2),S) = temp;
                
                temp = juicer2mat(...
                    juicer_tools_dump_mat('observed','KR',...
                    [selpath.hic,'/',samples.hic{S},'/inter_30.hic'],...
                    chr,chr,'BP',1E5));
                temp(isnan(temp)) = 0;
                H.(fields{i}).intra.juicer_KR{chr}(1:size(temp,1),1:size(temp,2),S) = temp;
            end
        end
    end
end
end

%% Eigenspectrum and VNGE - OLD
if 1==0
%%%%%%%%%%%%%%%%%
parameters.num_eigs = 10;
parameters.H_norm_method = 'corr';
%%%%%%%%%%%%%%%%%

fields = fieldnames(H);
for i = 1:numel(fields)
    fprintf('extracting eig-spectrum and VNGE Hi-C %s\n',fields{i})
    
    H.(fields{i}).vn_entropy = zeros(height(H.(fields{i}).chr_info), length(samples.hic));
    H.(fields{i}).eig_spec = zeros(height(H.(fields{i}).chr_info), parameters.num_eigs,length(samples.hic));
    
    %% calculate VNGE and eigenspectrum
    for chr = 1:22
        fprintf('extracting eig-spectrum and VNGE Hi-C %s chr%d\n',fields{i},chr)
        [H.(fields{i}).vn_entropy(chr,:), D_all] =...
            hicVnEntropy(double(H.(fields{i}).intra.juicer_oe_trim{chr}),...
            parameters.num_eigs, 1, parameters.H_norm_method);
        H.(fields{i}).eig_spec(chr,:,:) = reshape(D_all,1,parameters.num_eigs,length(samples.hic));
    end
    
    %% raw eigenspec of OE 100kb matrices
    if i==2
        D = zeros(parameters.num_eigs,22,length(samples.names));
        for chr = 1:22
            for S = 1:length(samples.names)
                D(:,chr,S) = eigs(corr(double(H.s100kb.intra.juicer_oe_trim{chr}(:,:,S))),...
                    parameters.num_eigs,'LM');
                D(:,chr,S) = D(:,chr,S)/sum(D(:,chr,S));
            end
        end
        for S = 1:length(samples.names)
            temp_color = hsv(10);
            
            figure('position',[100 100 800 400])
            b = bar(D(:,:,S)','stacked','FaceColor','flat');
            for eig_color = 1:parameters.num_eigs
                b(eig_color).CData = repmat(temp_color(eig_color,:),size(b(eig_color).CData,1),1);
            end
            title(sprintf('%s Normalized Eigenvalues',strrep(samples.names{S},'_',' ')))
            box off
            set(gca,'LineWidth',2)
            set(gca,'FontSize',15)
            ylim([0 1])
            saveas(gcf,sprintf('%s/%s_100kb_norm_%d_eig_barplot.png',...
                selpath.out,samples.names{S},parameters.num_eigs))
            % Eigenvalue decomposition of each 100kb intra-chromosomal
            % matrix, for all samples. Matrices are KR and O/E
            % normalized (See Methods), unmappable regions are removed, and
            % the eigenvalues of the correlation matrix are calculated. The
            % top 10 eigenvalues for each matrix are taken normalize to 1
            % in the plots
        end
    end
    
    %% VNE difference
    parameters.vne_max_diff = max(max(H.(fields{i}).vn_entropy,[],2)-min(H.(fields{i}).vn_entropy,[],2));
    for S1 = 1:length(samples.names)
        for S2 = S1+1:length(samples.names)
            figure('position',[100 100 800 400])
            bar(H.(fields{i}).vn_entropy(1:22,S1)-H.(fields{i}).vn_entropy(1:22,S2))
            box off
            set(gca,'LineWidth',2)
            set(gca,'FontSize',15)
            title(sprintf('%s - %s VNE Difference',strrep(samples.names{S1},'_',' '),...
                strrep(samples.names{S2},'_',' ')))
            ylabel('Entropy Difference')
            saveas(gcf,sprintf('%s/%s_v_%s_100kb_vneDiff_barplot.png',...
                selpath.out,samples.names{S1},samples.names{S2}))
            %ylim([-parameters.vne_max_diff parameters.vne_max_diff])
            
            % Entropy difference (VNE) between samples for each chromosome.
            % Matrices are KR and O/E normalized (See Methods), unmappable
            % regions are removed, and the eigenvalues of the correlation
            % matrix are calculated. The top 10 eigenvalues for each matrix
            % are normalized to 1. VNE is calculated from this
            % eigen-spectrum
        end
    end
    
    
    %% figure
    if 1==0
        %     figure, bar(1:22,H.(fields{i}).eig_spec(1:22,:,1),'stacked')
        %     title('normalized eigenspectrum, per chr,')
    end
end
end







% this script with calclulate the eigenspectrum and von neumann entropy of
% different cell types, and compare between them

clear
close all

restoredefaultpath
addpath(genpath('\\172.17.109.24\internal_4dn\tools\matlab\functions'))
addpath(genpath('\\172.17.109.24\internal_4dn\tools\matlab\scott_stephen_SF_functions'))

%% format data
CT = {'ALDH_P';'ALDH_N';'FIB';'ESC'};
hic_loc = {...
    '\\172.17.109.24\internal_4dn\projects\Wicha_BCSC\processed\hic\juicer\topdir_Sample_69655\inter_30.hic';...
    '\\172.17.109.24\internal_4dn\projects\Wicha_BCSC\processed\hic\juicer\topdir_Sample_69656\inter_30.hic';...
    '\\172.17.109.24\internal_4dn\projects\Fibroblast_ts\processed\hic\juicer_analysis_S34030\inter_30.hic';...
    '\\172.17.109.24\internal_4dn\projects\public_data\4DNFIOX3BGNE.hic'};
    
info = table(CT,hic_loc);

%% load Hi-C
%%% PARAMETERS vvv
hic_param.bin_type = 'BP';
hic_param.bin_size = 1E5;
hic_param.norm_1d = 'KR';
hic_param.norm_3d = 'oe';
%%% PARAMETERS ^^^

H = cell(22,1);
H2 = cell(22,1);
for Sample = 1:height(info)
    for Chr = 1:22
        fprintf('Sample:%s (%d/%d), chr:%d...\n',info.CT{Sample},Sample,height(info),Chr)
        if Sample==4
            temp = juicer2mat(juicer_tools_dump_mat(hic_param.norm_3d,hic_param.norm_1d,...
                info.hic_loc{Sample},Chr,Chr,hic_param.bin_type,hic_param.bin_size));
            H2{Chr} = padconcatenation_sr(H{Chr},temp,3); %concatenates matrices of different sizes
        else
            % .hic to matlab format
            temp = juicer2mat(juicer_tools_dump_mat(hic_param.norm_3d,hic_param.norm_1d,...
                info.hic_loc{Sample},Chr,Chr,hic_param.bin_type,hic_param.bin_size));
            H{Chr} = padconcatenation_sr(H{Chr},temp,3); %concatenates matrices of different sizes
        end
        
    end
end

%% trim Hi-C
%%% PARAMETERS vvv
hic_param.num_diag = 2;
hic_param.num_sparse = 0;
%%% PARAMETERS ^^^

H_trim = cell(22,1);
H2_trim = cell(22,1);
for Chr = 1:22
    fprintf('Trimming Hi-C Chr:%d...\n',Chr)
    H2_trim{Chr} = hic_trim(H2{Chr},hic_param.num_diag,hic_param.num_sparse);
    H_trim{Chr} = hic_trim(H{Chr},hic_param.num_diag,hic_param.num_sparse);
end

%% calculate eig spec and VNE
%%% PARAMETERS vvv
hic_param.num_eigs = 20;
hic_param.H_norm_method = 'corr';
%%% PARAMETERS ^^^

eig_spec = zeros(height(info),22,hic_param.num_eigs);
VNE = zeros(height(info),22);

for Sample = 1:height(info)
    for Chr = 1:22
        fprintf('Calculating eigspec and VNE. Sample:%s (%d/%d), chr:%d...\n',...
            info.CT{Sample},Sample,height(info),Chr)
        
        if Sample == 4
            % calculate eig spec
            [VN_entropy, D_all] = hic_VN_entropy(H2_trim{Chr}(:,:,Sample),...
                hic_param.num_eigs, 1, hic_param.H_norm_method);
        else
            % calculate eig spec
            [VN_entropy, D_all] = hic_VN_entropy(H_trim{Chr}(:,:,Sample),...
                hic_param.num_eigs, 1, hic_param.H_norm_method);
        end
        
        VNE(Sample,Chr) = VN_entropy;
        eig_spec(Sample,Chr,:) = D_all;
    end
end

%% figure eigspec
for Sample = 1:height(info)    
    [h] = plot_hic_eigspec(squeeze(eig_spec(Sample,:,:)),info.CT{Sample});
    saveas(gcf,sprintf('./figures/eigspec/eigspec_%s.png',info.CT{Sample}))
end

%% figure VNE compare
for Sample1 = 1:height(info)
    for Sample2 = Sample1+1:height(info)
        [h] = plot_hic_vne_diff([VNE(Sample1,:);VNE(Sample2,:)],...
            {info.CT{Sample1};info.CT{Sample2}});
        saveas(gcf,sprintf('./figures/vne/vne_diff_%s_%s.png',...
            info.CT{Sample1},info.CT{Sample2}))
    end
end

%% raw images
chr_select = 14;

figure
subplot(1,3,1), imagesc(log(H{chr_select}(:,:,1))), axis square
subplot(1,3,2), imagesc(log(H{chr_select}(:,:,3))), axis square
subplot(1,3,3), imagesc(log(H2{chr_select}(:,:,4))), axis square
linkaxes(get(gcf,'children'))


%% EXTRA (delete later)
% 
% temp = juicer2mat(juicer_tools_dump_mat(hic_param.norm_3d,hic_param.norm_1d,...
%     info.hic_loc{1},1,2,hic_param.bin_type,1E6));


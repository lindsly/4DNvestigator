function HiC_density_norm = norm_hic_bins(hic,bin_info,bin_name,Nbins,plot_options,mode,norm_unit)
%norm_hic_bins normalizes Hi-C matrices based on bin size and location
    %hic: NxNxT Hi-C matrices
    %bin_info: chromosome location of nodes. [chr,start,end]
    %bin_name: name of bins
    %Nbins: # of bins for histcounts
    %plot_options: optional plot if == 0,
    %mode: 
    %norm_unit: normaliization unit, default 1e5
    
    %example
    %load('\\172.17.109.24\internal_4DN\tools\matlab\data\TP0_HiC_genelevel.mat');
    %load('\\172.17.109.24\internal_4DN\tools\matlab\data\GeneTADinfo.mat', 'GeneInfo')
    %idx = 1:100;
    %hic = gene_contact(idx,idx);
    %gene_locs = cell2mat(GeneInfo(:,4));
    %gene_locs = gene_locs(idx,:);
    %A = norm_hic_bins(hic,gene_locs);

if nargin<7;norm_unit=1e5;end
if nargin<6;mode=0;end
if nargin<5;plot_options.fig=0;end
if nargin<4;Nbins=5e2;end
if nargin<3;bin_name=cell(size(hic,1),1);end

T = size(hic,3);
bin_info(:,2:3) = bin_info(:,2:3)/norm_unit;
Gene_len_sel = abs(bin_info(:,3)-bin_info(:,2));
Gene_cord_lin_sel = (bin_info(:,3)+bin_info(:,2))/2;
Gene_cord_chr_sel = bin_info(:,1);

Gene_len_mat = Gene_len_sel*Gene_len_sel.';
Gene_dist_chr_mat = abs( Gene_cord_chr_sel*ones(1,length(Gene_cord_chr_sel)) - ones(length(Gene_cord_chr_sel),1)*Gene_cord_chr_sel.' );
Gene_dist_chr_mat = ( Gene_dist_chr_mat < 1e-3 );
Gene_dist_lin_mat = abs( Gene_cord_lin_sel*ones(1,length(Gene_cord_lin_sel)) - ones(length(Gene_cord_lin_sel),1)*Gene_cord_lin_sel.' );
Gene_dist_lin_mat_tmp = Gene_dist_lin_mat;
Gene_dist_lin_mat_tmp(~Gene_dist_chr_mat) = nan;
Gene_dist_lin_mat_tmp(find(eye(length(Gene_dist_lin_mat_tmp)))) = nan;
Gene_dist_vec = Gene_dist_lin_mat(find(~tril(ones(size(Gene_dist_lin_mat))) & Gene_dist_chr_mat )); %%% upper triangluar part
[Countsbin,Bdbins] = histcounts(Gene_dist_vec,Nbins); %% similar distance within same chromosomes

HiC_density_norm = hic; %eps = 1e-10;
val_bins = zeros(Nbins,T);
text_legend = cell(T,1);
for t = 1:T
    text_legend{t} = sprintf('t%d',t);
    HiCtmp = hic(:,:,t);
    HiC_RPM_tmp = HiCtmp/1e3; %%% contacts at all of chrs per time; RPM
    HiC_density_tmp = HiC_RPM_tmp./Gene_len_mat; %%% contact frequency
    HiC_density_tmp = 0.5*(HiC_density_tmp+HiC_density_tmp.');
    HiC_density_norm_tmp = HiC_density_tmp;
    
    for i = 1:Nbins
        lb_tmp = Bdbins(i);
        up_tmp = Bdbins(i+1);
        if mode == 0
            mu_temp = mean ( nonzeros( HiC_density_tmp ( ( Gene_dist_lin_mat_tmp >= lb_tmp ) & ( Gene_dist_lin_mat_tmp < up_tmp ) ) ) ); %%% intra chrs
        else
            mu_temp = mean ( ( HiC_density_tmp ( ( Gene_dist_lin_mat_tmp >= lb_tmp ) & ( Gene_dist_lin_mat_tmp < up_tmp ) ) ) );
        end
        if isnan(mu_temp)
            mu_temp = 0; %%% no contact;
        end
        
        val_bins(i,t) = mu_temp;
        
        HiC_density_norm_tmp ( ( Gene_dist_lin_mat_tmp >= lb_tmp ) & ( Gene_dist_lin_mat_tmp < up_tmp ) ) = ...
            HiC_density_tmp ( ( Gene_dist_lin_mat_tmp >= lb_tmp ) & ( Gene_dist_lin_mat_tmp < up_tmp ) )/...
            mu_temp ; %%% intra-chr: return nan if no nonzero entries
    end
    
    if mode == 0
        tmp1 = mean( nonzeros(HiC_density_tmp( ( isnan(Gene_dist_lin_mat_tmp) & (~eye(length(Gene_dist_lin_mat_tmp))) ) )));
    else
        tmp1 = mean((HiC_density_tmp( ( isnan(Gene_dist_lin_mat_tmp) & (~eye(length(Gene_dist_lin_mat_tmp))) ) )));
    end
    
    HiC_density_norm_tmp(isnan(Gene_dist_lin_mat_tmp) & (~eye(length(Gene_dist_lin_mat_tmp))) ) ...
        = HiC_density_tmp(isnan(Gene_dist_lin_mat_tmp) & (~eye(length(Gene_dist_lin_mat_tmp))) )/...
        tmp1; %%% inter-chr part, return nan if no zero entries
    
    HiC_density_norm_tmp(isnan(HiC_density_norm_tmp))=0;
    
    if sum(sum( xor(HiC_density_norm_tmp ,HiCtmp ) )) > 0
        error('error in normalization');
    end
    
    HiC_density_norm(:,:,t) = HiC_density_norm_tmp;
    if plot_options.fig == 1
        cmap = [1,1,1; 1,.98,.98; 1,.96,.96; 1,.93,.93; 1,.9,.9; 1,.86,.86; 1,.8,.8;1,.6,.6; ...
            1,.4,.4; 1,.2,.2; 1,.1,.1; 1,.05,.05; 1,.02,.02; 1,0,0 ];
        h_tmp = figure;
        HiC_density_norm_tmp = (HiC_density_norm_tmp + HiC_density_norm_tmp.')*0.5;
        HiC_density_norm_tmp = HiC_density_norm_tmp - diag(diag(HiC_density_norm_tmp));
        imagesc(nthroot(log2(HiC_density_norm_tmp+0.5)+1,4));
        hold on; colormap(cmap) ;
        hcolor = [];
        hcolor(1) = plot(NaN, NaN, 'Marker','s','MarkerEdgeColor','r'...
            ,'MarkerFaceColor','r','MarkerSize',5,'LineStyle','none');
        axis tight;
        str_legend = [];
        str_legend{1} = sprintf('%d',full(ceil(log2((max(max(HiC_density_norm_tmp)))))));
        legend(hcolor,str_legend);
        title(sprintf('%s Norm Time%d',plot_options.title,t));
        savefig(h_tmp,fullfile(plot_options.folder, sprintf('%s_t%d',plot_options.title,t)));close(h_tmp);
    end
    fprintf('normalization t = %d\n',t);
end

end
function HiC_density_norm = normalization_GeneHiC_PM(HiC,Gene_names,GeneInfo,Nbins, plot_options,mode)
    
    T = size(HiC,3);
    idx_genes = cellfun(@(x) find(strcmp(x,GeneInfo(:,1))),Gene_names,'UniformOutput',false);
    idx_genes = cell2mat(idx_genes);
    Gene_Info_sel = GeneInfo(idx_genes,:);
    Gene_loc_sel = cell2mat(Gene_Info_sel(:,4));
    %%% 100Kb gene loc coordinate
    Gene_loc_sel(:,2:3) = Gene_loc_sel(:,2:3)/1e5; %%% 100Kb unit
    Gene_len_sel = abs(Gene_loc_sel(:,3)-Gene_loc_sel(:,2));
    Gene_cord_lin_sel = (Gene_loc_sel(:,3)+Gene_loc_sel(:,2))/2;
    Gene_cord_chr_sel = Gene_loc_sel(:,1);
    
    Gene_len_mat = Gene_len_sel*Gene_len_sel.';
    Gene_dist_chr_mat = abs( Gene_cord_chr_sel*ones(1,length(Gene_cord_chr_sel)) - ones(length(Gene_cord_chr_sel),1)*Gene_cord_chr_sel.' );
    Gene_dist_chr_mat = ( Gene_dist_chr_mat < 1e-3 );
    Gene_dist_lin_mat = abs( Gene_cord_lin_sel*ones(1,length(Gene_cord_lin_sel)) - ones(length(Gene_cord_lin_sel),1)*Gene_cord_lin_sel.' );
    Gene_dist_lin_mat_tmp = Gene_dist_lin_mat; 
    Gene_dist_lin_mat_tmp(~Gene_dist_chr_mat) = nan;
    Gene_dist_lin_mat_tmp(find(eye(length(Gene_dist_lin_mat_tmp)))) = nan;
    Gene_dist_vec = Gene_dist_lin_mat(find(~tril(ones(size(Gene_dist_lin_mat))) & Gene_dist_chr_mat )); %%% upper triangluar part
    [Countsbin,Bdbins] = histc(Gene_dist_vec,Nbins); %% histcounts similar distance within same chromosomes
    
    HiC_density_norm = HiC; %eps = 1e-10;
    val_bins = zeros(Nbins,T);
    text_legend = cell(T,1);
    for t = 1:T
        text_legend{t} = sprintf('t%d',t);
        HiCtmp = HiC(:,:,t);
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

%         if mode == 0
%            tmp = mean(nonzeros(diag(HiC_density_tmp)))*eye(length(HiC_density_tmp)) + ones(length(HiC_density_tmp)) - eye(length(HiC_density_tmp));  
%         else
%            tmp = mean((diag(HiC_density_tmp)))*eye(length(HiC_density_tmp)) + ones(length(HiC_density_tmp)) - eye(length(HiC_density_tmp));  
%         end
%         
%         HiC_density_norm_tmp = HiC_density_norm_tmp./tmp;

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
%             thr_cut = quantile(vec(HiC_density_norm_tmp(HiC_density_norm_tmp>0)),1-1e-3);
%             HiC_density_norm_tmp(HiC_density_norm_tmp>thr_cut) = thr_cut;
%             imagesc(nthroot(HiC_density_norm_tmp,4)); 
            HiC_density_norm_tmp = (HiC_density_norm_tmp + HiC_density_norm_tmp.')*0.5;
            HiC_density_norm_tmp = HiC_density_norm_tmp - diag(diag(HiC_density_norm_tmp));
            %imagesc(HiC_density_norm_tmp>0+0); 
            imagesc(nthroot(log2(HiC_density_norm_tmp+0.5)+1,4));
            hold on; colormap(cmap) ; 
             hcolor = [];
             hcolor(1) = plot(NaN, NaN, 'Marker','s','MarkerEdgeColor','r'...
                ,'MarkerFaceColor','r','MarkerSize',5,'LineStyle','none');
            axis tight; %colorbar;
            str_legend = [];
            str_legend{1} = sprintf('%d',full(ceil(log2((max(max(HiC_density_norm_tmp)))))));
            legend(hcolor,str_legend);
            title(sprintf('%s Norm Time%d',plot_options.title,t));
            savefig(h_tmp,fullfile(plot_options.folder, sprintf('%s_t%d',plot_options.title,t)));close(h_tmp);
%            figure,imagesc(HiC_density_norm_tmp>0), colormap hot, colormap(flipud(colormap));
        end
        disp(sprintf('normalization t = %d',t));
    end
%     if plot_options.fig == 1
%         h_tmp = figure;  %plot(val_bins);  legend(text_legend);
%         semilogy(0.5*(Bdbins(1:end-1)+Bdbins(2:end)),mean(val_bins,2)+1e-10,'LineStyle','-','Marker','.','MarkerSize',6,'Color','b','MarkerFaceColor','r','MarkerEdgeColor','r');
% %         shadedErrorBar([],mean(val_bins,2),0.8*std(val_bins,0,2),{'r-o','markerfacecolor','b'})
% %         set(gca,'YScale','log');
%         xlabel('Dist. between two genes (100 kb unit)');
%         ylabel('mean contact density');
%         title(sprintf('MeanContact versus dist: %s',plot_options.title));
%         savefig(h_tmp,fullfile(plot_options.folder, sprintf('MeanContact_versus_dist_%s',plot_options.title)));close(h_tmp);
%     end

end
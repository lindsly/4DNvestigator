function [] = networkExamples(Folder_Data, Folder_Result)
    %%% MyoD dataset
    load(fullfile(Folder_Data, 'GeneTADinfo.mat')); %%% gene info.
    load(fullfile(Folder_Data, 'MyoD_gene_Rna_raw.mat')); %%% RNAseq
    load(fullfile(Folder_Data, 'MyoD_Mb_HiC_raw.mat')); %%% HiC 1mb; wrong time order

    run_samps_old = {'64585','64584','71530','71533','71535','71537','71538','71536'...
        '71534','71532','71529','71531'};
    run_samps = {'64584','64585','71530','71529','71533','71532','71531','71536'...
        '71535','71534','71538','71537'};
    
    [LIA,time_order] = ismember(run_samps,run_samps_old) ;
    numReads = numReads(time_order);
    HiC_raw_1Mb = HiC_raw(:,:,time_order);
    T = min([size(HiC_raw,3),size(FPKM_gene,2)]); %%% T = 12
    RNAseq_AllT_AllRep = FPKM_gene_indiv(:,1:3*T); %%% replicates
    RNAseq_AllT = FPKM_gene(:,1:T); %%% ave. rnaseq
    Gene_Names_RNAseq = gene_names; %%% gene_names corresponding to RNAseq
    load(fullfile(Folder_Data, 'ChromosomeNumbers.mat')); %%% HiC 1mb
    Chr_start_MyoD_1Mb_raw = chr;

    %%% MYOD1 & MYOG
    genes_MYODMYOG = {'MYOD1','MYOG'};
    idx_genes = cell2mat(cellfun(@(x) find(strcmp(x,Gene_Names_RNAseq)),genes_MYODMYOG,'UniformOutput',false));
    RNAseq_MYODMYOG = RNAseq_AllT(idx_genes,:);
    t_real = [-48,0:8:80];
    Tvec = 2:12; %%% 3: 8hr ; 4-16, 5-24, 6-32, 7-40,8-48,9-56
    for t = 1:length(Tvec)
        t_real_str{t} = sprintf('%d',t_real(Tvec(t)));
    end

    %%%%% color library for plotting
    color_set = [rgb('Black');rgb('DarkSlateGray');rgb('Red');rgb('Salmon');rgb('Crimson');...
        rgb('DarkRed');rgb('Pink');rgb('HotPink');rgb('DeepPink');rgb('Orange');...
        rgb('Tomato');rgb('OrangeRed');rgb('Gold');rgb('DarkKhaki');rgb('Tan');rgb('RosyBrown');...
        rgb('Chocolate');rgb('Green');rgb('LightGreen');rgb('YellowGreen');rgb('Olive');rgb('Teal');...
        rgb('Blue');rgb('Cyan');rgb('DeepSkyBlue');rgb('Purple');rgb('BlueViolet');rgb('Indigo')];

    cmap = [1,1,1; 1,.98,.98; 1,.96,.96; 1,.93,.93; 1,.9,.9; 1,.86,.86; 1,.8,.8;1,.6,.6; ...
        1,.4,.4; 1,.2,.2; 1,.1,.1; 1,.05,.05; 1,.02,.02; 1,0,0 ];

    color_set_time = [[1 0 0];[0.87,0.49,0];[0.25,0.25,0.25];[0.75,0.75,0];[1, 0.6,0.78];[0 1 0];...
        [0 1 1]; [0,0.45,0.74];[0 0 1]; [0.49,0.18,0.56];[1 0 1];[0.64,0.08,0.18]];

    %%% multilayer network: gene network, two groups

    %% multilayer network
    %%% Target genes: skin + muscle (myotube)
    load(fullfile(Folder_Data, sprintf('go_terms_formatted.mat')));
    Gene_Names_skin_muscle = unique(union(go.fibroblast,go.myotube));
    Gene_Names_skin_muscle = unique(union(Gene_Names_skin_muscle,go.myoblast));

    %%% Fib Gnetwork
    load(fullfile(Folder_Data, 'MyoD_GHiC_raw_all_Fib.mat'));
    GeneInfo_sort_Fib = GeneInfo_sort;
    idx_sel_Fib_tmp = cell2mat( cellfun(@(x) find(strcmp(x,GeneInfo_sort_Fib(:,1))), Gene_Names_skin_muscle,'UniformOutput',false) );
    idx_sel_Fib_tmp = sort(idx_sel_Fib_tmp,'ascend');
    GeneInfo_Fib = GeneInfo_sort_Fib(idx_sel_Fib_tmp,:);
    HiC_G_Fib = HiC_G_Raw( idx_sel_Fib_tmp,idx_sel_Fib_tmp,: );

    %%% MyoD Gnetwork
    load(fullfile(Folder_Data, 'MyoD_GHiC_raw_skin_muscle.mat'));
    GeneInfo_sort_MyoD = GeneInfo_sort;
    idx_sel_MyoD_tmp = cell2mat( cellfun(@(x) find(strcmp(x,GeneInfo_sort_MyoD(:,1))), Gene_Names_skin_muscle,'UniformOutput',false) );
    idx_sel_MyoD_tmp  = sort(idx_sel_MyoD_tmp,'ascend');
    GeneInfo_MyoD = GeneInfo_sort_MyoD(idx_sel_MyoD_tmp,:);
    HiC_G_MyoD = HiC_G_Raw( idx_sel_MyoD_tmp,idx_sel_MyoD_tmp,: );

    %%% common gene set
    Gene_common =  intersect(GeneInfo_Fib(:,1),GeneInfo_MyoD(:,1)); %%% sort;
    idx_sel_tmp1 = cell2mat( cellfun(@(x) find(strcmp(x,GeneInfo_Fib(:,1))), Gene_common,'UniformOutput',false) );
    idx_sel_tmp1 = sort(idx_sel_tmp1,'ascend');
    GeneInfo_Fib_common = GeneInfo_Fib(idx_sel_tmp1,:);
    HiC_G_Fib_common = HiC_G_Fib(idx_sel_tmp1,idx_sel_tmp1,:);
    coord_chr_Fib = []; loc_gene_temp = cell2mat(GeneInfo_Fib_common(:,4));
    for i_chr = 1:22
        idx_temp = find(loc_gene_temp(:,1) == i_chr);
        if ~isempty(idx_temp)
            coord_chr_Fib = [coord_chr_Fib; [idx_temp(1),idx_temp(end)] ];
        else
            coord_chr_Fib = [coord_chr_Fib; [NaN,NaN] ];
        end
    end

    idx_sel_tmp2 = cell2mat( cellfun(@(x) find(strcmp(x,GeneInfo_MyoD(:,1))), Gene_common,'UniformOutput',false) );
    idx_sel_tmp2 = sort(idx_sel_tmp2,'ascend');
    GeneInfo_MyoD_common = GeneInfo_MyoD(idx_sel_tmp2,:);
    HiC_G_MyoD_common = HiC_G_MyoD(idx_sel_tmp2,idx_sel_tmp2,:);
    coord_chr_MyoD = []; loc_gene_temp = cell2mat(GeneInfo_MyoD_common(:,4));
    for i_chr = 1:22
        idx_temp = find(loc_gene_temp(:,1) == i_chr);
        if ~isempty(idx_temp)
            coord_chr_MyoD = [coord_chr_MyoD; [idx_temp(1),idx_temp(end)] ];
        else
            coord_chr_MyoD = [coord_chr_MyoD; [NaN,NaN] ];
        end
    end

    for j = 1:size(GeneInfo_Fib_common,1)
        if ~strcmp(GeneInfo_Fib_common{j,1}, GeneInfo_MyoD_common{j,1})
            error('error in name mismatch');
        end
    end
    GeneInfo_Fib_MyoD_common = GeneInfo_MyoD_common;

    %%% Gene expression
    %%% Fib
    load(fullfile(Folder_Data, 'TsFib_kb_HiC_GeneTADinfo.mat'));
    RNAseq_AllT_Fib = FIBGeneRNAseq(:,1:8);
    Gene_Names_RNAseq_Fib = GeneInfo(:,1);
    idx_sel_tmp1 = cell2mat( cellfun(@(x) find(strcmp(x,Gene_Names_RNAseq_Fib)), GeneInfo_Fib_MyoD_common(:,1),'UniformOutput',false) );
    RNAseq_Fib = RNAseq_AllT_Fib(idx_sel_tmp1,:);

    %%% MyoD
    idx_sel_tmp2 = cell2mat( cellfun(@(x) find(strcmp(x,Gene_Names_RNAseq)), GeneInfo_Fib_MyoD_common(:,1),'UniformOutput',false) );
    RNAseq_MyoD = RNAseq_AllT(idx_sel_tmp2,:);

    %%% HiC normalization
    Nbins = 5e2; mode = 0;
    plot_options.fig = 0;
    plot_options.title = sprintf('Fib G2G HiC norm');
    plot_options.folder = Folder_Result;

    thr_min_norm = 0.1;thr_max_norm = 1e-2;

    %% normalization for Fib
    HiC_G_Fib_common_norm = normalization_GeneHiC_PM(HiC_G_Fib_common,GeneInfo_Fib_common(:,1),GeneInfo_Fib_common,Nbins, plot_options,mode); %% RPM
    for t = 1:size(HiC_G_Fib_common_norm,3)
        Ht = HiC_G_Fib_common_norm(:,:,t);
        Ht = Ht - diag(diag(Ht));
        Ht = 0.5*(Ht + Ht.');

        thr_cut_max = quantile(vec(Ht(Ht>0)),1-thr_max_norm);
        thr_cut_min = quantile(vec(Ht(Ht>0)),thr_min_norm);
        Ht(Ht>Ht) = thr_cut_max;
        Ht(Ht<thr_cut_min) = 0;

        HiC_G_Fib_common_norm(:,:,t) =Ht;
    end
    %% normalization for MyoD
    plot_options.title = sprintf('MyoD G2G HiC norm');
    HiC_G_MyoD_common_norm = normalization_GeneHiC_PM(HiC_G_MyoD_common,GeneInfo_MyoD_common(:,1),GeneInfo_MyoD_common,Nbins, plot_options,mode); %% RPM
    for t = 1:size(HiC_G_MyoD_common_norm,3)
        Ht = HiC_G_MyoD_common_norm(:,:,t);
        Ht = Ht - diag(diag(Ht));
        Ht = 0.5*(Ht + Ht.');

        %%%
        thr_cut_max = quantile(vec(Ht(Ht>0)),1-thr_max_norm);
        thr_cut_min = quantile(vec(Ht(Ht>0)),thr_min_norm);
        Ht(Ht>Ht) = thr_cut_max;
        Ht(Ht<thr_cut_min) = 0;
        HiC_G_MyoD_common_norm(:,:,t) =Ht;
    end

    %%% multiplex network: Fib & MyoD seperately, write in one loop
    t_real_all = [ 0    8    16    24    32    40    48    56];
    Group_t = [0 8 16 24; 32 40 48 56];
    dataset_name = {'MYOD','Fib'};

    deg_alpha_binary = cell(size(Group_t,1),length(dataset_name));
    odeg_binary = cell(size(Group_t,1),length(dataset_name));
    p_muli_binary= cell(size(Group_t,1),length(dataset_name));
    eps = 1e-2;
    for i_dataset = 1:length(dataset_name)
        for i_group = 1:size(Group_t,1)
            t_real = Group_t(i_group,:);
            [~,idx_t_real,~] = intersect(t_real_all,t_real);
            if i_dataset == 1    %%% MyoD
                Tvec = 2:9; %%% 0hr.. 2:9
                Tvec_MyoD = Tvec(idx_t_real);
                T = length(Tvec_MyoD);
                HiC_multiplex = HiC_G_MyoD_common_norm(:,:,Tvec_MyoD);
            else %%% Fib
                Tvec = 1:8; %%% 0hr.. 2:9
                Tvec_Fib = Tvec(idx_t_real);
                T = length(Tvec_Fib);
                HiC_multiplex = HiC_G_Fib_common_norm(:,:,Tvec_Fib);
            end
            for alpha = 1:size(HiC_multiplex,3)
                H_alpha = HiC_multiplex(:,:,alpha);
                H_alpha = 0.5*(H_alpha + H_alpha.');
                H_alpha = H_alpha - diag(diag(H_alpha));
                H_alpha_binary = (H_alpha > 1e-5) + 0;
                deg_alpha_binary{i_group,i_dataset}(:,alpha) = sum(H_alpha_binary,2) + eps;
            end %%% multilayer
            odeg_binary{i_group,i_dataset} = sum(deg_alpha_binary{i_group,i_dataset},2);
            %%% multiplex coefficient
            temp_data_binary = deg_alpha_binary{i_group,i_dataset}./((odeg_binary{i_group,i_dataset})*ones(1,size(deg_alpha_binary{i_group,i_dataset},2)));
            p_muli_binary{i_group,i_dataset} = (1 - sum(temp_data_binary.^2,2))*T/(T-1);
        end
    end

    %%% post processing; focused on myod
    idx_sel = find( ( sum(odeg_binary{1,1},2) + sum(odeg_binary{2,1},2) ) > 5  );
    GeneInfo_Fib_MyoD_sel = GeneInfo_Fib_MyoD_common(idx_sel,:);
    thr_sel_var = 2;
    odeg_binary_norm = cell(size(Group_t,1),length(dataset_name));
    p_muli_binary_norm = cell(size(Group_t,1),length(dataset_name));
    for i_dataset = 1:length(dataset_name)

        odeg_binary_sel_temp = [odeg_binary{1,i_dataset}(idx_sel),odeg_binary{2,i_dataset}(idx_sel)];
        p_muli_binary_sel_temp = [p_muli_binary{1,i_dataset}(idx_sel),p_muli_binary{2,i_dataset}(idx_sel)];
        odeg_binary_sel_temp_norm = ( odeg_binary_sel_temp-min(odeg_binary_sel_temp(:)) )./( max(odeg_binary_sel_temp(:))-min(odeg_binary_sel_temp(:)) );
        p_muli_binary_sel_temp_norm = ( p_muli_binary_sel_temp-min(p_muli_binary_sel_temp(:)) )./( max(p_muli_binary_sel_temp(:))-min(p_muli_binary_sel_temp(:)) );
        odeg_binary_norm{1,i_dataset} = odeg_binary_sel_temp_norm(:,1);
        odeg_binary_norm{2,i_dataset} = odeg_binary_sel_temp_norm(:,2);
        p_muli_binary_norm{1,i_dataset} = p_muli_binary_sel_temp_norm(:,1);
        p_muli_binary_norm{2,i_dataset} = p_muli_binary_sel_temp_norm(:,2);
    end

    %%%% MYOD data
    diff_str = sqrt( ( odeg_binary_norm{1,1} - odeg_binary_norm{2,1} ).^2 + ( p_muli_binary_norm{1,1} - p_muli_binary_norm{2,1} ).^2 );
    thr_str = mean(diff_str)+thr_sel_var*std(diff_str);
    idx_sel_str = find(diff_str>=thr_str);
    [~,idx_sort]= sort(diff_str(idx_sel_str),'descend');
    idx_sel_str = idx_sel_str(idx_sort);
    hfig = figure('name','Degree vs Multiplex Participation Coefficient - All Genes');
    plot(p_muli_binary_norm{1,1},odeg_binary_norm{1,1} ,'Marker','^', 'MarkerSize',4,'MarkerFaceColor','m','MarkerEdgeColor','m','LineStyle','none'); hold on;
    plot(p_muli_binary_norm{2,1},odeg_binary_norm{2,1}  ,'Marker','s', 'MarkerSize',4,'MarkerFaceColor','g','MarkerEdgeColor','g','LineStyle','none'); hold on;
    plot(p_muli_binary_norm{1,1}(idx_sel_str) ,odeg_binary_norm{1,1}(idx_sel_str),'Marker','^','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','none','LineStyle','none'); hold on;
    plot(p_muli_binary_norm{2,1}(idx_sel_str),odeg_binary_norm{2,1}(idx_sel_str) ,'Marker','s','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','none','LineStyle','none'); hold on;
    hcolor = [];   %
    hcolor(1) = plot(NaN,NaN,'Marker','^', 'MarkerSize',8,'MarkerFaceColor','m','MarkerEdgeColor','m','LineStyle','none'); hold on;
    hcolor(2) = plot(NaN,NaN,'Marker','s', 'MarkerSize',8,'MarkerFaceColor','g','MarkerEdgeColor','g','LineStyle','none'); hold on;
    hcolor(3) = plot(NaN,NaN,'Marker','o','MarkerSize',9,'MarkerEdgeColor','k','MarkerFaceColor','none','LineStyle','none'); hold on;
    hcolor(4) = plot(NaN,NaN,'Marker','s','MarkerSize',9,'MarkerEdgeColor','k','MarkerFaceColor','none','LineStyle','none'); hold on;
    legend(hcolor,{'0-24 hrs', '32-56 hrs','Selected genes (0-24 hrs)','Selected genes (32-56 hrs)'});  legend('boxoff');
    ax = gca;
    ax.XDir = 'reverse';
    xlabel('Multiplex participation coeff.');
    ylabel('Overlapping degree');

    %%%% shift
    hfig = figure('name','Degree vs Multiplex Participation Coefficient - Selected Genes');
    plot(p_muli_binary_norm{1,1}(idx_sel_str) ,odeg_binary_norm{1,1}(idx_sel_str) ,'Marker','^', 'MarkerSize',4,'MarkerFaceColor','m','MarkerEdgeColor','m','LineStyle','none'); hold on;
    plot(p_muli_binary_norm{2,1}(idx_sel_str),odeg_binary_norm{2,1}(idx_sel_str),'Marker','s', 'MarkerSize',4,'MarkerFaceColor','g','MarkerEdgeColor','g','LineStyle','none'); hold on;
    plot([p_muli_binary_norm{1,1}(idx_sel_str).';p_muli_binary_norm{2,1}(idx_sel_str).'],[odeg_binary_norm{1,1}(idx_sel_str).';odeg_binary_norm{2,1}(idx_sel_str).'],'LineStyle','--','Color',rgb('darkgrey'),'LineWidth',1); hold on;
    text(p_muli_binary_norm{1,1}(idx_sel_str)-0.03,odeg_binary_norm{1,1}(idx_sel_str),GeneInfo_Fib_MyoD_sel(idx_sel_str,1)); hold on;
    hcolor = [];   %
    hcolor(1) = plot(NaN,NaN,'Marker','^', 'MarkerSize',8,'MarkerFaceColor','m','MarkerEdgeColor','m','LineStyle','none'); hold on;
    hcolor(2) = plot(NaN,NaN,'Marker','s', 'MarkerSize',8,'MarkerFaceColor','g','MarkerEdgeColor','g','LineStyle','none'); hold on;
    hcolor(3) = plot(NaN,NaN,'LineStyle','--','Color',rgb('darkgrey'));
    legend(hcolor,{'Selected genes (0-24 hrs)','Selected genes (32-56 hrs)','Gene pairs'});  legend('boxoff');
    ax = gca;
    ax.XDir = 'reverse';
    xlabel('Multiplex participation coefficient');
    ylabel('Overlapping degree');

    %%% Network point of view: interpret this phenomenan
    Genes_star = GeneInfo_Fib_MyoD_sel;
    idx_genes_temp = cell2mat(cellfun(@(x) find(strcmp(x,GeneInfo_Fib_MyoD_common(:,1))),Genes_star(:,1),'UniformOutput',false));
    Genes_network = GeneInfo_Fib_MyoD_common(idx_genes_temp,:);

    HiC_multiplex = zeros(length(idx_genes_temp),length(idx_genes_temp),2*T); iter = 1;
    idx_MyoD = []; idx_Fib = [];

    for t = 1:(size(Group_t,1)*size(Group_t,2))
        HiC_multiplex(:,:,iter) = HiC_G_MyoD_common_norm(idx_genes_temp,idx_genes_temp,t);
        iter = iter + 1;
    end
    HiC_multiplex_cum = sum(HiC_multiplex,3);
    HiC_multiplex_cum = HiC_multiplex_cum - diag(diag(HiC_multiplex_cum));
    HiC_multiplex_cum = 0.5*(HiC_multiplex_cum+HiC_multiplex_cum.');
    HiC_multiplex_cum_bi = (HiC_multiplex_cum > 0) + 0;

    %%% selected genes
    GeneInfo_Fib_MyoD_common_deg_nz_sel=GeneInfo_Fib_MyoD_sel(idx_sel_str(1:3),:);
    idx_genes_temp2 = cell2mat(cellfun(@(x) find(strcmp(x,Genes_network(:,1))),GeneInfo_Fib_MyoD_common_deg_nz_sel(:,1),'UniformOutput',false));

    HiC_multiplex_sel_cum_bi_sel = HiC_multiplex_cum_bi( idx_genes_temp2,:);
    idx_gene_clean = find( sum(HiC_multiplex_sel_cum_bi_sel,1) > 1e-3 ); %%% rectangular matrix
    idx_gene_clean = unique(union(idx_gene_clean,idx_genes_temp2));

    Genes_network_clean = Genes_network(idx_gene_clean,:);
    HiC_multiplex_clean = HiC_multiplex(idx_gene_clean,idx_gene_clean,:);
    HiC_multiplex_cum_bi_clean = HiC_multiplex_cum_bi(idx_gene_clean,idx_gene_clean);

    BGobj = biograph(HiC_multiplex_cum_bi_clean, Genes_network_clean(:,1), 'ShowArrows','off','ShowWeights','off',...
        'EdgeType','straight','LayoutType','equilibrium');
    dolayout(BGobj);

    %%% extract 2D location
    xy_loc = [];
    for i_node = 1:size(HiC_multiplex_cum_bi_clean,1)
        xy_loc = [xy_loc;BGobj.nodes(i_node).Position];
    end
    sz_Genes_tar_temp = cell2mat(Genes_network_clean(:,4));
    sz_Genes_tar = abs(sz_Genes_tar_temp(:,3)-sz_Genes_tar_temp(:,2))+1;
    [~, idx_node_sort]= sort(sz_Genes_tar,'descend');

    %%% size of node
    num_node_sz_level = 5;
    node_weight = sz_Genes_tar./max(sz_Genes_tar);
    node_weight_quantile = quantile(unique(node_weight),linspace(0.1, 1, num_node_sz_level));

    %%% edges
    color_muliplex = {'r','b','g','m'};
    edge_linspec = linspace(0.5, 0.5, num_node_sz_level);
    node_sz_spec = linspace(5, 5, num_node_sz_level);
    color_edge = [ [255 210 210]/255;[210 210 255]/255 ;[210 255 210]/255 ;[255 204 249]/255];
    color_node = [153 153 153]/255;
    color_text = [0 0 0]/255;

    idx_genes_sig_temp = cell2mat(cellfun(@(x) find(strcmp(x,Genes_network_clean(:,1))),GeneInfo_Fib_MyoD_common_deg_nz_sel(:,1),'UniformOutput',false));
    GroupTnew = [1,2,3,4;5,6,7,8];
    for j = 1:size(GroupTnew,1)
        T_group = GroupTnew(j,:); lenT = length(T_group);
        hfig = figure('name',['Network Visualization ', num2str(j)]);
        for i = 1:lenT
            subplot(2,2,i);
            Adj_t = HiC_multiplex_clean(:,:,T_group(i));
            Adj_t = Adj_t - diag(diag(Adj_t));
            Adj_t = 0.5*(Adj_t + Adj_t.');

            H_triu = triu(Adj_t,1);
            [ki,kj]=find(H_triu);
            edge_weight = H_triu(find(H_triu));
            edge_weight = edge_weight./max(edge_weight);
            edge_weight_quantile = quantile(edge_weight,linspace(0.1, 1, num_node_sz_level));

            for ii = 1:length(edge_weight_quantile)
                if ii == 1
                    ind_edge = find( edge_weight  <= edge_weight_quantile(ii) );
                else
                    ind_edge = find( edge_weight  <= edge_weight_quantile(ii) & edge_weight  > edge_weight_quantile(ii-1) );
                end
                ki_temp = ki(ind_edge);
                kj_temp = kj(ind_edge);
                Lia = find(ismember(ki_temp,idx_genes_sig_temp));
                Lib = find(ismember(kj_temp,idx_genes_sig_temp));
                plot([xy_loc(ki_temp,1),xy_loc(kj_temp,1)].',...
                    [xy_loc(ki_temp,2),xy_loc(kj_temp,2)].',...
                    '-', 'LineWidth',edge_linspec(ii),...
                    'Color', color_edge(i,:)); hold on;

                if ~isempty(Lia)
                    ki_temp2 = ki_temp(Lia);
                    kj_temp2 = kj_temp(Lia);
                    plot([xy_loc(ki_temp2,1),xy_loc(kj_temp2,1)].',...
                        [xy_loc(ki_temp2,2),xy_loc(kj_temp2,2)].',...
                        '-', 'LineWidth',edge_linspec(ii),...
                        'Color', color_muliplex{i}); hold on;
                end

                if ~isempty(Lib)
                    ki_temp3 = ki_temp(Lib);
                    kj_temp3 = kj_temp(Lib);
                    plot([xy_loc(ki_temp3,1),xy_loc(kj_temp3,1)].',...
                        [xy_loc(ki_temp3,2),xy_loc(kj_temp3,2)].',...
                        '-', 'LineWidth',edge_linspec(ii),...
                        'Color', color_muliplex{i});  hold on;
                end

            end %% edge weights


            for i_node = 1:length(idx_node_sort)

                idx_node = idx_node_sort(i_node);
                xy_loc_nodei = xy_loc(idx_node,:); %%% xy loc
                weight_nodei = node_weight(idx_node);

                idx_weight = find( weight_nodei  >=  node_weight_quantile);
                if isempty(idx_weight)
                    sz_nodei = node_sz_spec(1);
                else
                    sz_nodei = node_sz_spec(idx_weight(end));
                end

                if  ismember(idx_node,idx_genes_sig_temp)
                    plot(xy_loc_nodei(1), xy_loc_nodei(2), 'Marker', 'o', 'MarkerEdgeColor', 'm', ...
                        'MarkerFaceColor', 'm','MarkerSize',sz_nodei); hold on;
                    text(xy_loc_nodei(1), xy_loc_nodei(2),Genes_network_clean{idx_node,1},'Color',color_text,'FontWeight','bold'); hold on;
                else
                    plot(xy_loc_nodei(1), xy_loc_nodei(2), 'Marker', 'o', 'MarkerEdgeColor', color_node, ...
                        'MarkerFaceColor', color_node,'MarkerSize',sz_nodei); hold on;
                end

            end %%% nodes
            hold on;   axis tight;  set(gca, 'XTick', []); set(gca, 'YTick', []);
            title(sprintf('%d hr',t_real_all(T_group(i))));
        end
    end


    %%% inter-chr or intra-chr
    %%% Genes_network_clean

    %%% selected genes
    GeneInfo_Fib_MyoD_common_deg_nz_sel=GeneInfo_Fib_MyoD_sel(idx_sel_str,:);
    idx_genes_temp2 = cell2mat(cellfun(@(x) find(strcmp(x,Genes_network(:,1))),GeneInfo_Fib_MyoD_common_deg_nz_sel(:,1),'UniformOutput',false));

    HiC_multiplex_sel_cum_bi_sel = HiC_multiplex_cum_bi( idx_genes_temp2,:);
    idx_gene_clean = find( sum(HiC_multiplex_sel_cum_bi_sel,1) > 1e-3 ); %%% rectangular matrix
    idx_gene_clean = unique(union(idx_gene_clean,idx_genes_temp2));

    Genes_network_clean = Genes_network(idx_gene_clean,:);
    HiC_multiplex_clean = HiC_multiplex(idx_gene_clean,idx_gene_clean,:);
    HiC_multiplex_cum_bi_clean = HiC_multiplex_cum_bi(idx_gene_clean,idx_gene_clean);

    Deg_in_out_chrs = zeros(size(Genes_network_clean,1),3,2);  %%% in out total for PM
    ratio_MFdiff_in_out_chrs = zeros(size(Genes_network_clean,1),2); %%% in & out PM diff
    num_MFdiff_in_out_chrs = zeros(size(Genes_network_clean,1),2);
    loc_gene_temp = cell2mat(Genes_network_clean(:,4));
    eps0 = 1e-5;  coord_chr = [];
    for i_chr = 1:23
        idx_temp = find(loc_gene_temp == i_chr);
        if ~isempty(idx_temp)
            coord_chr = [coord_chr; [idx_temp(1),idx_temp(end)] ];
        else
            coord_chr = [coord_chr; [NaN,NaN] ];
        end
    end

    for ig = 1:size(Genes_network_clean,1)
        Genei = Genes_network_clean(ig,1);
        %%% find genei's location
        idx_Genei = cell2mat(cellfun(@(x) find(strcmp(x,Genes_network_clean(:,1))),Genei,'UniformOutput',false));
        chr_gi = loc_gene_temp(idx_Genei,1);
        idx_gene_inChr = [coord_chr(chr_gi,1):coord_chr(chr_gi,2)];
        if isnan(idx_gene_inChr)
            idx_gene_inChr = [];
        end
        idx_gene_outChr = setdiff( (1:size(Genes_network_clean,1)), idx_gene_inChr);
        %%%G1
        Edge_multiplex_M_inChr = sum(vec(HiC_multiplex_clean(idx_Genei,idx_gene_inChr,GroupTnew(1,:))>eps0))/(length(idx_gene_inChr));
        Edge_multiplex_M_outChr = sum(vec(HiC_multiplex_clean(idx_Genei,idx_gene_outChr,GroupTnew(1,:))>eps0))/(length(idx_gene_outChr));
        Edge_multiplex_M_total = sum(vec(HiC_multiplex_clean(idx_Genei,:,GroupTnew(1,:))>eps0))/(length(idx_gene_inChr)+length(idx_gene_outChr));

        Deg_in_out_chrs(ig,:,1) = [Edge_multiplex_M_inChr,Edge_multiplex_M_outChr,Edge_multiplex_M_total];

        %%%G2
        Edge_multiplex_F_inChr = sum(vec(HiC_multiplex_clean(idx_Genei,idx_gene_inChr,GroupTnew(2,:))>eps0))/(length(idx_gene_inChr));
        Edge_multiplex_F_outChr = sum(vec(HiC_multiplex_clean(idx_Genei,idx_gene_outChr,GroupTnew(2,:))>eps0))/(length(idx_gene_outChr));
        Edge_multiplex_F_total = sum(vec(HiC_multiplex_clean(idx_Genei,:,GroupTnew(2,:))>eps0))/(length(idx_gene_inChr)+length(idx_gene_outChr));

        Deg_in_out_chrs(ig,:,2) = [Edge_multiplex_F_inChr,Edge_multiplex_F_outChr,Edge_multiplex_F_total];

        %%% PM in/out Chr
        ratio_MFdiff_in_out_chrs(ig,1) = ...
            (Edge_multiplex_M_inChr-Edge_multiplex_F_inChr) / (Edge_multiplex_M_total-Edge_multiplex_F_total);
        ratio_MFdiff_in_out_chrs(ig,2) = ...
            (Edge_multiplex_M_outChr-Edge_multiplex_F_outChr) / (Edge_multiplex_M_total-Edge_multiplex_F_total);

        num_MFdiff_in_out_chrs(ig,1) = (Edge_multiplex_M_inChr-Edge_multiplex_F_inChr);
        num_MFdiff_in_out_chrs(ig,2) = (Edge_multiplex_M_outChr-Edge_multiplex_F_outChr);
    end

    ratio_MFdiff_in_out_chrs(find(isnan(ratio_MFdiff_in_out_chrs))) = 0;
    ratio_MFdiff_in_out_chrs(find(isinf(ratio_MFdiff_in_out_chrs))) = 0;

    %%% Fig: Plot InChr outChr of selected genes
    genes_specific = GeneInfo_Fib_MyoD_sel(idx_sel_str,1);
    idx_genes_sig_temp = cell2mat(cellfun(@(x) find(strcmp(x,Genes_network_clean(:,1))),genes_specific,'UniformOutput',false));

    idx_sel = idx_genes_sig_temp;
    gene_name_sel =[];
    for ii = 1:length(idx_sel)
        temp = cell2mat(Genes_network_clean(idx_sel(ii),4));
        gene_name_sel{ii} = sprintf('%d: %s',temp(1),Genes_network_clean{idx_sel(ii),1});
    end

    hfig = figure('name','Inter- and Intra-Chromosomal Contacts');

    plot(1:length(idx_sel), Deg_in_out_chrs(idx_sel,1,1),'Marker','o','MarkerSize',7, 'LineStyle','-','Color','r','MarkerFaceColor','r','MarkerEdgeColor','r' );
    hold on;
    plot(1:length(idx_sel), Deg_in_out_chrs(idx_sel,1,2),'Marker','o','MarkerSize',7, 'LineStyle','-','Color','r','MarkerFaceColor','none','MarkerEdgeColor','r' );
    hold on;
    plot(1:length(idx_sel), Deg_in_out_chrs(idx_sel,2,1),'Marker','s','MarkerSize',7, 'LineStyle','--','Color','b','MarkerFaceColor','b','MarkerEdgeColor','b' );
    hold on;
    plot(1:length(idx_sel), Deg_in_out_chrs(idx_sel,2,2),'Marker','s','MarkerSize',7, 'LineStyle','--','Color','b','MarkerFaceColor','none','MarkerEdgeColor','b' );
    hold on;     hcolor = [];      str_legend = [];

    hcolor(1) = plot(NaN,NaN,'Marker','o','markersize',7,'LineStyle','-','Color','r','MarkerFaceColor','r','MarkerEdgeColor','r'); hold on;
    hcolor(2) = plot(NaN,NaN,'Marker','o','markersize',7,'LineStyle','-','Color','r','MarkerFaceColor','none','MarkerEdgeColor','r'); hold on;
    hcolor(3) = plot(NaN,NaN,'Marker','s','markersize',7,'LineStyle','--','Color','b','MarkerFaceColor','b','MarkerEdgeColor','b'); hold on;
    hcolor(4) = plot(NaN,NaN,'Marker','s','markersize',7,'LineStyle','--','Color','b','MarkerFaceColor','none','MarkerEdgeColor','b'); hold on;

    str_legend{1} = sprintf('Intra-chr contacts, Reprogramming');
    str_legend{2} = sprintf('Intra-chr contacts, Proliferation');
    str_legend{3} = sprintf('Inter-chr contacts, Reprogramming');
    str_legend{4} = sprintf('Inter-chr contacts, Proliferation');

    axis tight;
    set(gca,'xtick',1:length(idx_sel)) ;
    set(gca,'xticklabel',gene_name_sel);   xtickangle(90);
    legend(hcolor,str_legend);
    ylabel('Contact density per gene');


    %% Save figures
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
      FigHandle = FigList(iFig);
      FigName   = get(FigHandle, 'Name');
      savefig(FigHandle, [Folder_Result, '\',FigName, '.fig']);
      saveas(FigHandle, [Folder_Result, '\',FigName, '.png']);
    end
end

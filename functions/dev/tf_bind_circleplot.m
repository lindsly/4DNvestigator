function [h] = tf_bind_circleplot(TF_names,method,std_amount,bin_size)
%tf_bind_circleplot plots the location of TF
%   Detailed explanation goes here

if nargin<4;bin_size = 1E6;end
if nargin<3;std_amount = 1;end
if nargin<2;method = 'DGC';end
if nargin<1;TF_names = {'SOX2','POU5F1','MYC','KLF4','NFKB1'};end

switch method
    case 'DGC'
        load('\\172.17.109.24\internal_4DN\tools\matlab\data\GeneTFBinding_547_v2.mat',...
            'TFConsensusMotifs','TFNames')      % TF names infor
        load('\\172.17.109.24\internal_4DN\tools\matlab\data\GeneTADinfo.mat',...
            'GeneInfo')                         % Gene loc info
        load('\\172.17.109.24\internal_4DN\tools\matlab\data\hg_chr_lengths.mat')   % chr lengths
        
        Ydat = cell(22,length(TF_names));
        Gene_locs = cell2mat(GeneInfo(:,4));
        idx = find_ismember_ordered(TFNames,TF_names);
        
        %% gene TSS to chr bin
        gene_tss_1mb = zeros(length(GeneInfo),2);
        for i = 1:length(GeneInfo)
            if strcmp(GeneInfo{i,5},'+')
                gene_tss_1mb(i,:) = [Gene_locs(i,1),ceil(Gene_locs(i,2)/bin_size)];
            else
                gene_tss_1mb(i,:) = [Gene_locs(i,1),ceil(Gene_locs(i,3)/bin_size)];
            end
        end
        
        %% create table
        info = table(GeneInfo(:,1),gene_tss_1mb,TFConsensusMotifs(:,idx));
        info = splitvars(info);
        info.Properties.VariableNames = [{'gene_name','chr','mb_bin'},TF_names];
        
        %% filter, sort by loc  
        thresh = repmat(mean(info{:,4:end})+std_amount*std(info{:,4:end}),height(info),1);
        info{:,4:end} = info{:,4:end}>thresh;
        info = sortrows(info,{'chr','mb_bin'});
        
        chr_length = ceil(chr_length.grch37(1:22)/bin_size);
        
        for chr_i = 1:length(chr_length)            
            temp = info(info.chr==chr_i,3:end);
            temp2_stats = grpstats(temp,'mb_bin','max');
            for i = 1:length(TF_names)
                Ydat{chr_i,i} = zeros(chr_length(chr_i),1);
                Ydat{chr_i,i}(temp2_stats.mb_bin) = temp2_stats{:,2+i};
                Ydat{chr_i,i} = logical(Ydat{chr_i,i});
            end
        end
        
        %% TF_names locs
        idx2 = find_ismember_ordered(GeneInfo(:,1),TF_names);
        TF_loc_info = [Gene_locs(idx2,1),ceil(Gene_locs(idx2,2:3)/bin_size)];
        
end

%% Geoff Prettiness!
%%%%%%%%%%%%%%%%%%%%%%
labelspace = 0;
r       = 0.9;
R       = 1;
pad     = chr_length(1)*.1;%20;
% J       = data2color(1:22,[softblue;softgreen;softyellow;softred]);
arroww  = 2.5/360*2*pi;
color_scheme = distinguishable_colors(length(TF_names)).*0.8+.2;
%%%%%%%%%%%%%%%%%%%%%%

cs = chr_length;%cs = chrmsizes; cs=cs(:,4);
N  = sum(cs)+(22-sign(labelspace))*pad+sign(labelspace);

chrm_start = cumsum(cs)-cs + pad*(0:21)' + 1;
chrm_end   = chrm_start + cs;

theta   = pi/2+linspace(2*pi*(1-labelspace),2*pi*labelspace,N);
x       = r*cos(theta);
y       = r*sin(theta);
X       = R*cos(theta);
Y       = R*sin(theta);
h       = cell(22,1);
clf;

if 1==1
    temp = linspace(0.55, 0.85, length(TF_names)+1);
    yamr = temp(2:end);
    yamR = temp(1:end-1); temptheta=linspace(0,2*pi,800);
else
    yamr = 0.55:0.1:0.85;
    yamR = 0.45:0.1:0.75; temptheta=linspace(0,2*pi,800);
end
for i=1:length(yamr);
    plot((yamr(i)+yamR(i))/2*cos(temptheta),(yamr(i)+yamR(i))/2*sin(temptheta),'-','color',0.8*[1 1 1]); hold on;
end
plot(0.95*cos(temptheta),0.95*sin(temptheta),'-','color',0.8*[1 1 1]);
for c=1:22
    h{c}=patch([x(chrm_start(c):chrm_end(c)) X(chrm_end(c):-1:chrm_start(c))],...
        [y(chrm_start(c):chrm_end(c)) Y(chrm_end(c):-1:chrm_start(c))],0.9*[1 1 1],'EdgeColor',0.5*[1 1 1],'LineWidth',2); hold on;
    set(h{c},'ButtonDownFcn',@(src,~) geoffyo(c,src));
    if 1==0
    if c==3;
        plotData(ones(5,1),theta(chrm_start(c)+163),theta(chrm_start(c)+197),0.9,1,softred,0,1);         % 181
        plotDataArrow('SOX2',theta(chrm_start(c)+181)-arroww,theta(chrm_start(c)+181)+arroww,1.01,1.1,softred);
    elseif c==6;
        plotData(ones(5,1),theta(chrm_start(c)+14),theta(chrm_start(c)+48),0.9,1,softyellow,0,1);        % 31
        plotDataArrow('OCT4',theta(chrm_start(c)+31)-arroww,theta(chrm_start(c)+31)+arroww,1.01,1.1,softyellow);
    elseif c==8;
        plotData(ones(5,1),theta(chrm_start(c)+110),theta(chrm_start(c)+144),0.9,1,softgreen,0,1);       % 127
        plotDataArrow('MYC',theta(chrm_start(c)+127)-arroww,theta(chrm_start(c)+127)+arroww,1.01,1.1,softgreen);
    elseif c==9;
        plotData(ones(5,1),theta(chrm_start(c)+90),theta(chrm_start(c)+124),0.9,1,softblue,0,1);        % 107
        plotDataArrow('KLF4',theta(chrm_start(c)+107)-arroww,theta(chrm_start(c)+107)+arroww,1.01,1.1,softblue);
    end
    else
        if ismember(c,TF_loc_info(:,1))
            temp = find(ismember(TF_loc_info(:,1),c));
            temp2 = chrm_end(c)-chrm_start(c);
            for i = 1:length(temp)
                temp_s = TF_loc_info(temp(i),2)-1;
                temp_e = TF_loc_info(temp(i),3)+1;
                temp_m = round(mean(TF_loc_info(temp(i),2:3),2));
                plotData(ones(5,1),theta(chrm_start(c)+temp_s),theta(chrm_start(c)+temp_e),...
                    0.9,1,color_scheme(temp(i),:),0,1);         % 181
                plotDataArrow(TF_names{temp(i)},theta(chrm_start(c)+temp_m)-arroww,...
                    theta(chrm_start(c)+temp_m)+arroww,1.01,1.1,color_scheme(temp(i),:));
            end
        end
    end
    middle = round(mean([chrm_start(c);chrm_end(c)]));
    text(mean([x(middle);X(middle)]),...
        mean([y(middle);Y(middle)]),...
        num2str(c),'FontSize',8,'horizontalalignment','center','verticalalignment','middle','Rotation',...
        (-pi/2+mean([theta(chrm_start(c));theta(chrm_end(c))]))/2/pi*360);
end

for c=1:22;
    for i = 1:length(TF_names)
        plotData(1-Ydat{c,i},theta(chrm_start(c)),theta(chrm_end(c)),yamr(i),yamR(i),brighten(color_scheme(i,:),0.95),0,1);
        plotData(Ydat{c,i},theta(chrm_start(c)),theta(chrm_end(c)),yamr(i),yamR(i),color_scheme(i,:),0,1);
    end
    htemp=patch([x(chrm_start(c):chrm_end(c))/r*0.45 X(chrm_end(c):-1:chrm_start(c))/R*0.85],...
        [y(chrm_start(c):chrm_end(c))/r*0.45 Y(chrm_end(c):-1:chrm_start(c))/R*0.85],[1 1 1],'EdgeColor',0.5*[1 1 1],'LineWidth',1.25);
    set(htemp,'facecolor','none');
end

axis equal
axis tight
axis fill
axis off

hold off;
set(gcf,'color','w');

text(0,0.95,'Chromosome','fontsize',11,'horizontalalignment','center','verticalalignment','middle','backgroundcolor','w');
for i = 1:length(TF_names)
    text(0,(yamR(i)+yamr(i))/2,TF_names{i},'fontsize',11,'horizontalalignment',...
        'center','verticalalignment','middle','backgroundcolor','w');
    %text(0,(yamR(i)+yamr(i))/2,TF_names{i},'fontsize',11,'horizontalalignment',...
    %'center','verticalalignment','middle','backgroundcolor','none');
end
text(0,0,TF_names,'fontsize',25,'horizontalalignment','center','verticalalignment','middle');

end


%% subfunctions
function plotData(data,theta1,theta2,r1,r2,colors,fmin,fmax,fzero,outlineon)
if nargin<6 || isempty(colors);
    colors = sweetblue;
end
if nargin<7 || isempty(fmin);
    fmin = min(data);
end
if nargin<8 || isempty(fmax);
    fmax = max(data);
end
if nargin<9 || isempty(fzero);
    fzero = 0;
end
if nargin<10 || isempty(outlineon);
    outlineon=0;
end

zeroline = (fzero-fmin)/(fmax-fmin);
data  = (data'-fmin)./(fmax-fmin);

n=length(data);
idx   = [(1:n);(1:n)]; idx=idx(:);
data  = data(idx);
tops  = max(data,zeroline);
bots  = min(data,zeroline);

idx   = [(1:n+1);(1:n+1)]; idx=idx(2:end-1);
theta = linspace(theta1,theta2,n+1);
theta = theta(idx);

xd    = (((r2-r1)*bots)+r1).*cos(theta);
Xd    = (((r2-r1)*tops)+r1).*cos(theta);
yd    = (((r2-r1)*bots)+r1).*sin(theta);
Yd    = (((r2-r1)*tops)+r1).*sin(theta);

if outlineon
    patch([xd Xd(end:-1:1)],[yd Yd(end:-1:1)],[1 1 1],'edgecolor',[1 1 1],'linewidth',5);
end
if size(colors,1)==1
    patch([xd Xd(end:-1:1)],[yd Yd(end:-1:1)],colors,'edgecolor','none');
else
    for i=1:n
        patch([xd([2*i-1,2*i]) Xd([2*i,2*i-1])],...
            [yd([2*i-1,2*i]) Yd([2*i,2*i-1])],...
            colors(i,:),'edgecolor','none');
    end
end
end

function plotDataArrow(textlabel,theta1,theta2,r1,r2,color)
theta = linspace(theta1,theta2,89);
r     = linspace(r1,r2,4);
tidx  = [45 89 60 60 30 30 1];
ridx  = [ 1  2  2  4  4  2 2];
patch(r(ridx).*cos(theta(tidx)),r(ridx).*sin(theta(tidx)),color,'edgecolor','none');
text((r2+0.01)*cos(theta(45)),(r2+0.01)*sin(theta(45)),textlabel,...
    'horizontalalignment','left','verticalalignment','middle','rotation',theta(45)/2/pi*360)
end

function geoffyo(c,src,~)
J = data2color(1:22,[softblue;softgreen;softyellow;softred]);
src.FaceColor=J(c,:);
end


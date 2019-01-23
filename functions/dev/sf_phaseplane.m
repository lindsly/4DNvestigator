function [h] = sf_phaseplane(R,H,bin_names,sample_names,num_connections,rna_seq_norm)
%sf_phaseplane plots the phase plane for structure and function
%   R is an n x m matrix with rna-seq values. n number of genes/bins, m
%   number of samples
%   H is an n x m matrix with FDV values. n number of genes/bins, m
%   number of samples
%   rna_seq_normis an option to normalize the RNA-seq values

if nargin < 6;rna_seq_norm='none';end
if nargin < 5;num_connections=10;end
if nargin < 4;sample_names=cell(size(R,2),1);end
if nargin < 3;bin_names=cell(size(R,1),1);end

%% normalize data if necessary
switch rna_seq_norm
    case 'none'
    case 'log2'
        R=log(R+.5)+1;
    case 'sqrt'
        R=sqrt(R);
end

%% determine distance
dist_all = zeros(size(R,1)*nchoosek(size(R,2),2),4);
count = 1;

for i = 1:size(R,2)
    for ii = i+1:size(R,2)
        for iii=1:size(R,1)
            dist_all(count,:) = [pdist([R(iii,i),H(iii,i);R(iii,ii),H(iii,ii)]),...
                i,ii,iii];
            count = count+1;
        end
    end
end

[~,I] = sort(dist_all(:,1),'descend');

%% figure
figure, hold on
for i = 1:size(R,2)
    plot(R(:,i),H(:,i),'*')
end

for i = 1:num_connections
    loc = dist_all(I(i),2:end);
    plot([R(loc(3),loc(1)),R(loc(3),loc(2))],...
        [H(loc(3),loc(1)),H(loc(3),loc(2))],'k-',...
        'linewidth',5)
    if ~isempty(bin_names{loc(3)})
        text(mean([R(loc(3),loc(1)),R(loc(3),loc(2))]),...
            mean([H(loc(3),loc(1)),H(loc(3),loc(2))]),...
            bin_names{loc(3)},'backgroundcolor','w')
    end
end
xlabel('RNA-seq')
ylabel('FDV')
if ~isempty(sample_names{1})
    legend(sample_names)
end

%% extra
% dist_mat = squareform(pdist(centroids));
% plot(R,H,'*'),lengend(sample_names)
% xlabel('RNA-seq')
% ylabel('FDV')

%% sample
if 1==0
    sf_phaseplane(rand(100,2),rand(100,2))
end

end


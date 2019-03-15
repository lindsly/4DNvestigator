function TAD_interval = TAD_Laplace_Sijia(H,sig0, ms0)
% TADs extraction via Laplacian segmentation

% From Laura's 'Extract_TADs_from_JIE.m':
%chr_choice = 14;        %chr selection
%TP = 1;                 %TP selection
%FN_thresh = 0.8;        %FN threshold (smaller = larger TADS, start at 0.8)
%ms0 = 3;                %smallest allowable TAD size
%MODE = 0;               %0 = start at O/E. 1 = start at whole chr
%eig_num = 2;            %2 = FN


%% Default parameters
if nargin<3; ms0=3; end
if nargin<2; sig0=.7; end

L = size(H,1);
Pos = [1, L];

%% Recursive splitting

spa = zeros(L,1);
spa(Pos)=1;

for i = 1 : length(Pos)-1
    
    % Block range
    idx = Pos(i):Pos(i+1)-1;
    
    if length(idx)>ms0
        % Sub-matrix
        SH = H(idx,idx);
        
        % Fiedler number and vector
%         [Fdv,Fdvl]=Fdvectors(SH);
        [Fdvl,Fdv] = Fiedlervec(SH);
        
        % If the Fiedler number of the block is small
        if Fdvl <= sig0
            % Continue to split
            sp = SubSplit(SH,sig0,ms0);
            
            % Mark boundaries
            spa(Pos(i)+find(sp>0)-1)=1;
        end
        
    end
    
end


posn = find(spa>0);

posn = MergeSmall_Sijia(posn,H,ms0);

TAD_interval = posn; 
TAD_interval = unique([TAD_interval; L+1]);

%     %%%%%DELETEv%%%%%%%%%062816_sr
%     Hn=ToepNorm(H0);
%     figure, imagesc(Hn,[0 10]), colormap gray
%     for i = 1:length(posn)-1
%         hold on
%         plot([posn(i),posn(i)],[posn(i),posn(i+1)],'r-')
%         plot([posn(i),posn(i+1)],[posn(i),posn(i)],'r-')
%         plot([posn(i),posn(i+1)],[posn(i+1),posn(i+1)],'r-')
%         plot([posn(i+1),posn(i+1)],[posn(i),posn(i+1)],'r-')
%     end
%     length(posn)
%     %%%%%DELETE^%%%%%%%%%

% figure,imagesc(H+dH,[0,6]),colormap gray, colormap(flipud(colormap))
% hold on
% for i = 2 : length(posn)
%     %    plot([1,size(HHnm,1)],[posn(i),posn(i)]-0.5,'--','linewidth',2)
%     plot([posn(i-1),posn(i)],[posn(i-1),posn(i-1)],'r-','linewidth',2)
%     plot([posn(i-1),posn(i)],[posn(i),posn(i)],'r-','linewidth',2)
%     plot([posn(i-1),posn(i-1)],[posn(i-1),posn(i)],'r-','linewidth',2)
%     plot([posn(i),posn(i)],[posn(i-1),posn(i)],'r-','linewidth',2)
%     
% end
% axis equal
% xlim([1250 1350])
% ylim([1250 1350])
% hold off
% figure, bar(FDV_i), xlim([229 491]), ylim([-.06 .06])

end


%% ============= Sub-functions =============

%% Sub-function 1: Recusive split
function sp = SubSplit(Ho,sig0, ms0)
% Recursively splitting a connection matrix via Fiedler value and vector



% Fiedler vector
% [Fdv,Fdvl]=Fdvectors(Ho);
[Fdvl,Fdv] = Fiedlervec(Ho);

if abs(Fdvl) < 1e-4 %%% disconnected; cannot trust Fdv
    [Acc_tmp,pcc_tmp] = largest_component_simple(Ho); 
    pcc_tmp = pcc_tmp + 0;
    pcc_tmp(find(pcc_tmp==0)) = -1;
    Fdv = pcc_tmp; %%% split based on connected components
    disp('Disconnected subgraph in TAD');
end

% Position of sign change (sub-block starting)
Pos = [1;find(sign(Fdv(2:end))-sign(Fdv(1:end-1))~=0)+1;size(Ho,1)];
Pos = unique(Pos); 
sp = zeros(size(Ho,1),1);
sp(Pos) = 1; %%% Starting point of TAD
% If Fiedler value is high enough
if  Fdvl > sig0+1e-5        %   +1e-5 for numerical stability
    sp = 1;
    return;
end

% For each sub-block
for i = 1 :  length(Pos)-1
    % Range
    idx = Pos(i):Pos(i+1)-1;

    % minimum sub-block size
    if length(idx)>ms0
       % Continue to split
       sp1 = SubSplit(Ho(idx,idx),sig0,ms0);
       % Mark bock boundary
       sp(Pos(i)+find(sp1>0)-1)=1;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
    end
end


end

%% Fiedler vector calculation
function [Fdv,Fdvl] = Fdvectors(H,NN)

H = (H +H')/2;
N = size(H,1);

if nargin < 2
    NN = 0;
end
dgr = sum(H);
dgr(dgr==0)=1;
DN = diag(1./sqrt(dgr));
L = DN*(diag(dgr)-H)*DN;

if NN == 1
    L = diag(dgr)-H;
end

L = (L+L')/2;
[V,D] = eigs(L,2,'SA');  % Or we can just use svd / eig
Fdv = V(:,end);
Fdvl =(D(end,end));
% 
% L = (L+L')/2;
% LI = eye(size(L))-L;
% [U1,D1]=eigs(LI,1,'LA');
% LI = LI - U1*U1'*2;
% LI = (LI+LI')/2;
% [Fdv,Fdvl] = eigs(LI,1,'LA');
% Fdvl = 1 - Fdvl;


if NN == 1;
    Fdvl = Fdvl/size(L,1)^0.3;
end
end

%% Merging small regions
function Pos=MergeSmall(posn,H)
Pos = posn;
Posr = Pos;
Find region only with 1 bin size
idx1 = find(Pos(2:end)-Pos(1:end-1)==1);
for i = 1 : length(idx1)
    cond1 = idx1(i)+1 <= length(Pos);
    cond2 = idx1(i)-1 >= 1;
    % Check this single bin this more similar to upstream or downstream bins.
    if idx1(i)+1 <= length(Pos)
% %         vrp = mean(H(Pos(idx1(i)),Pos(idx1(i))+1:Pos(idx1(i)+1))); 
        if idx1(i)+2 > length(Pos)
           vrp = mean(H(Pos(idx1(i)),Pos(idx1(i))+1 :  Pos(idx1(i)+1)   ));  
        else
           vrp = mean(H(Pos(idx1(i)),Pos(idx1(i))+1 :  ( Pos(idx1(i)+2)-1 )   )); 
        end
    end
    if idx1(i)-1 >= 1
        vrm = mean(H(Pos(idx1(i)),Pos(idx1(i)-1):Pos(idx1(i))-1));
    end
    if cond1&cond2  & vrm > vrp
        Posr(idx1(i)) = -100;
    elseif  cond1&cond2  & vrm < vrp
        Posr(idx1(i)+1) = -100;
    end
end
Pos(Posr==-100)=[];

%     Pos_tmp = unique( [Posr; length(H)+1] );
%     Pos_tmp_modify = Pos_tmp;
%     
%     for i = 1: ( length(Pos_tmp) - 1 )
%         pos_TAD_tmp = Pos_tmp(i): ( Pos_tmp(i+1) - 1 );
%         if i > 1 && i < ( length(Pos_tmp) - 1 ) %%% do merge
%             pos_TAD_tmp_above = Pos_tmp(i+1): ( Pos_tmp(i+2) - 1 );
%             pos_TAD_tmp_below = Pos_tmp(i-1): ( Pos_tmp(i) - 1 );
%             vr_above = mean(mean(H(pos_TAD_tmp,pos_TAD_tmp_above),2));
%             vr_below = mean(mean(H(pos_TAD_tmp,pos_TAD_tmp_below),2));
%             if vr_above > vr_below
%                 Pos_tmp_modify(i+1) = -100;
%             else
%                 Pos_tmp_modify(i) = -100;
%             end
%         end
%     end
%     Pos_tmp(Pos_tmp_modify==-100)=[];
%     Pos = Pos_tmp;



end

function Pos = MergeSmall_Sijia(posn,H,ms0)

Pos = posn;
for len_sm = 1: ms0
    Posr = Pos;
% %     find small bin of size len_sm
    idx1 = find(Pos(2:end)-Pos(1:end-1)==len_sm);
    for i = 1 : length(idx1)
        cond1 = idx1(i)+1 <= length(Pos);
        cond2 = idx1(i)-1 >= 1;
        if idx1(i)+1 <= length(Pos)
            if idx1(i)+2 > length(Pos)
               vrp = mean(mean(H( Pos(idx1(i)) : Pos(idx1(i)+1)-1,  Pos(idx1(i)+1)  : length(Pos)   )));  
            else
               vrp = mean(mean(H( Pos(idx1(i)) : Pos(idx1(i)+1)-1 , Pos(idx1(i)+1) :  ( Pos(idx1(i)+2)-1 )   ))); 
            end
        end
        if idx1(i)-1 >= 1
            vrm = mean(mean(H(  Pos(idx1(i)) : Pos(idx1(i)+1)-1 , Pos(idx1(i)-1):Pos(idx1(i))-1 ) ));
        end
        if cond1&cond2  & vrm > vrp
            Posr(idx1(i)) = -100;
        elseif  cond1&cond2  & vrm < vrp
            Posr(idx1(i)+1) = -100;
        end
    end
    Pos(Posr==-100)=[]; %%% update Pos

end

end
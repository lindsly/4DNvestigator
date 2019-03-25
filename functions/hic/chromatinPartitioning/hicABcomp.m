function [ABcomp,groupIdx,FdNum] = hicABcomp(A,method,rnaSeq,rnaSeqNorm,...
    chrDivide,plotFlag,fdvSplitMethod,chrDivideLoc)
%hicABcomp determines Hi-C A/B compartments
%   hicABcomp determines Hi-C A/B compartment partitioning through either the
%   Fiedler vector or the first principal component of the Hi-C correlation
%   matrix.
%
%   Inputs
%   A:              Hi-C matrix, typically normalized (NxN double; default:
%                   N/A)
%   method:         method for AB compartment partitioning (string;
%                   default: 'fiedler')
%   rnaSeq:         RNA-seq values associated with Hi-C matrix (Nx1
%                   double; default: zeros)
%   rnaSeqNorm:     RNA-seq normalization method (string; default: 'log2')
%   chrDivide:      divide chromsome if arms are incorectly partitioned
%                   (string; default: 'yes')
%   plotFlag:       plot chr once completed compartmentalization (string;
%                   default: 'yes')
%   chrDivideLoc:   predetermined location of chr centromere, relative to
%                   Hi-C matrix
%
%   Outputs
%   ABcomp:         AB designation for each compartment (Nx1 double)
%   groupIdx:       Index for partition groups (Nx1 double)
%   FdNum:          Fiedler number
%
%
%   Version 1.1 (3/8/19)
%   Written by: Scott Ronquist
%   Contact: scotronq@umich.edu
%   Contributors:
%   Created: 1/22/19
%   Revision History:
%   v1.0 (1/22/19)
%   * created script
%   v1.1 (3/8/19)
%   * added input argument for chrDivideLoc

%% default parameters
if ~exist('fdvSplitMethod','var')||isempty(fdvSplitMethod);fdvSplitMethod='sign';end
if ~exist('plotFlag','var')||isempty(plotFlag);plotFlag=0;end
if ~exist('chrDivide','var')||isempty(chrDivide);chrDivide='no';end
if ~exist('rnaSeqNorm','var')||isempty(rnaSeqNorm);rnaSeqNorm='log2';end
if ~exist('rnaSeq','var')||isempty(rnaSeq);rnaSeq=zeros(length(A),1);end
if ~exist('method','var')||isempty(method);method='fiedler';end
if ~exist('chrDivideLoc','var')||isempty(chrDivideLoc);chrDivideLoc=[];end

method = lower(method);
chrDivide = lower(chrDivide);
rnaSeqNorm = lower(rnaSeqNorm);
FdNum = [];

%% AB compartment method
switch method
    case 'fiedler'
        [~,Fdv,~,FdNum] = hicLaplacianFdv(A);
        ABcomp = Fdv;
    case 'pc1'
        temp = corr(A);
        temp(isnan(temp))=0;
        temp = temp-diag(diag(temp))+diag(ones(length(temp),1));
        [PC1,~] = eigs(temp,1);
        ABcomp = PC1;
end

%% split in 2 if separated by chr arm
switch chrDivide
    case 'no'
    case 'yes'
        
        if ~isempty(chrDivideLoc)
            fprintf('splitting AB analysis by chr arm, user input...\n')
            splitLoc = chrDivideLoc;
            [abComp1] = hicABcomp(A(1:splitLoc-1,1:splitLoc-1),...
                method,rnaSeq(1:splitLoc-1),rnaSeqNorm,'no');
            [abComp2] = hicABcomp(A(splitLoc:end,splitLoc:end),...
                method,rnaSeq(splitLoc:end),rnaSeqNorm,'no');
            ABcomp = [abComp1;abComp2];
            
        else
            % t-test if to determine if partitioning is strictly by chr arm
            abDiff = nan(length(ABcomp),1);
            h = nan(length(ABcomp),1);
            p = nan(length(ABcomp),1);
            for i = round(length(ABcomp)*.25):round(length(ABcomp)*.75)
                abDiff(i) = abs(mean(ABcomp(1:i-1))-mean(ABcomp(i:end)));
                [h(i),p(i)] = ttest2(ABcomp(1:i-1),ABcomp(i:end));
            end
            
            % set threshold for t-test
            switch method
                case 'fiedler'
                    thresh = 1E-20;
                case 'pc1'
                    thresh = 0;
            end
            
            % repartition each arm if incorrect
            if nanmin(p) <= thresh
                [~,splitLoc] = max(abDiff);
                fprintf('splitting AB analysis by chr arm, computer estimation...\n')
                [abComp1] = hicABcomp(A(1:splitLoc-1,1:splitLoc-1),...
                    method,rnaSeq(1:splitLoc-1),rnaSeqNorm,'no');
                [abComp2] = hicABcomp(A(splitLoc:end,splitLoc:end),...
                    method,rnaSeq(splitLoc:end),rnaSeqNorm,'no');
                ABcomp = [abComp1;abComp2];
            end
        end
end

%% rna-seq normalization
switch rnaSeqNorm
    case 'none'
    case 'log2'
        rnaSeq = log2(rnaSeq+1);
end

%% correlate with RNA-seq, if available
if mean(rnaSeq(ABcomp > 0)) < mean(rnaSeq(ABcomp < 0))
    ABcomp = -ABcomp;
end

%% determine partition
groupIdx = zeros(1,length(ABcomp));
switch method
    case 'pc1'
        groupIdx(ABcomp > 0) = 2;
        groupIdx(ABcomp < 0) = 1;
    case 'fiedler'
        switch fdvSplitMethod
            case 'bisection'
                fdvSplitVal = median(ABcomp);
                groupIdx(ABcomp > fdvSplitVal) = 2;
                groupIdx(ABcomp < fdvSplitVal) = 1;
            case 'ratio'
                %{
                https://digitalcommons.csbsju.edu/cgi/viewcontent.cgi?referer
                =https://www.google.com/&httpsredir=1&article=1021&context=honors_theses
                %}
                [ABcompSorted,sortIdx] = sort(ABcomp);
                ASorted = A(sortIdx,sortIdx);
                [A_,~] = ratioCut(ASorted);
                fdvSplitVal = mean(ABcompSorted(length(A_):length(A_)+1));
                groupIdx(ABcomp > fdvSplitVal) = 2;
                groupIdx(ABcomp < fdvSplitVal) = 1;
                
            case 'sign'
                fdvSplitVal = 0;
                groupIdx(ABcomp > fdvSplitVal) = 2;
                groupIdx(ABcomp < fdvSplitVal) = 1;
                
            case 'gap'
                ABcompSorted = sort(ABcomp);
                [~,maxIdx] = max(diff(ABcompSorted));
                fdvSplitVal = mean(ABcompSorted(maxIdx:maxIdx+1));
                groupIdx(ABcomp > fdvSplitVal) = 2;
                groupIdx(ABcomp < fdvSplitVal) = 1;
        end
end

%% plot (commented)
if plotFlag==1
    % RNA-seq, AB compartment, Hi-C, column sum
    figure, subplot(7,1,1),plot(rnaSeq),axis tight, ylabel(sprintf('RNAseq\n norm: %s',rnaSeqNorm))
    subplot(7,1,2),plot(ABcomp),axis tight, ylabel(sprintf('A/B compartment\n method: %s',method))
    subplot(7,1,3:6),plot_hic(log(A),'erez',[prctile(A(:),1) prctile(A(:),99)]), ylabel('Hi-C matrix')
    subplot(7,1,7),plot(sum(A)),axis tight, ylabel('Column sum')
    
    linkaxes(get(gcf,'children'),'x')
    
    % AB compartment bargraph with color
    figure, b = bar(ABcomp,'FaceColor','flat','EdgeColor','none');
    b.CData(groupIdx==1,:) = repmat([1 0 0],sum(groupIdx==1),1);
    b.CData(groupIdx==2,:) = repmat([0 1 0],sum(groupIdx==2),1);
end

end


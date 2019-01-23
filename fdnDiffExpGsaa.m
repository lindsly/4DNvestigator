function [R] = fdnDiffExpGsaa(dataInfo,R)
%fdnDiffExpGsaa This function performs differential expression and
%creates files for GSAA
%
%   Input
%   dataInfo: 4DNvestigator structure that contains information on samples input
%   R: 4DNvestigator structure that contains RNA-seq data
%
%   Output 
%   R: 4DNvestigator structure that contains RNA-seq data
%
%   Scott Ronquist, scotronq@umich.edu. 1/22/19

%% start
% create cell array with differential expression comparisons
samps = unique(dataInfo.sampleInfo.sample);
Tps = unique(dataInfo.sampleInfo.timePoint);
diffExpTable = cell(length(samps)*length(Tps));
nameAll = cell(length(samps)*length(Tps),1);

% create comparison matrix
[A,B] = meshgrid(1:length(samps),1:length(Tps));
c=cat(2,A',B');
compIndex = reshape(c,[],2);

% get output folder
fnBase = [dataInfo.path.output,dataInfo.delim,'data',dataInfo.delim,'gsaa',...
    dataInfo.delim];

% compare between all possibilities
for i = 1:size(compIndex,1)
    for ii = i+1:size(compIndex,1)
        
        % convert to real sample names
        tempSampleA = samps{compIndex(i,1)};
        tempTpA = Tps(compIndex(i,2));
        tempSampleB = samps{compIndex(ii,1)};
        tempTpB = Tps(compIndex(ii,2));
        
        % get A and B locs relative to R.TPM variable
        tempALocs = ismember(R.TPM.Properties.VariableNames,...
            dataInfo.sampleInfo.name(not(cellfun('isempty',...
            regexp(dataInfo.sampleInfo.uniqueName,...
            sprintf('rnaseq_s%s_t%i*',tempSampleA,tempTpA))))));
        tempBLocs = ismember(R.TPM.Properties.VariableNames,...
            dataInfo.sampleInfo.name(not(cellfun('isempty',...
            regexp(dataInfo.sampleInfo.uniqueName,...
            sprintf('rnaseq_s%s_t%i*',tempSampleB,tempTpB))))));
        
        % create geneTable
        meanA = mean(R.TPM{:,tempALocs},2);
        meanB = mean(R.TPM{:,tempBLocs},2);
        meanBase = (meanA + meanB) / 2;
        foldChange = meanA ./ meanB;
        log2FC = log2(foldChange);
        
        geneTable = table(meanBase,meanA,meanB,foldChange,log2FC);
        geneTable = [R.TPM(:,1:6),geneTable];
        
        [pvalue,padj] = matlabNegbinDE(R.expected_count{:,tempALocs}, R.expected_count{:,tempBLocs});
        geneTable.pvalue = pvalue;
        geneTable.padj = padj;
        
        % add to full comparison table
        diffExpTable{i,ii} = geneTable;
        
        %% GSAASeqSP analysis
        % http://gsaa.unc.edu/userguide_gsaaseqsp.html
        % general info
        FileName = sprintf('%s_vs_%s',...
            sprintf('rnaseq_s%s_t%i',tempSampleA,tempTpA),...
            sprintf('rnaseq_s%s_t%i',tempSampleB,tempTpB));
        FileNameFull = sprintf('%s%s',fnBase,FileName);
        
        % Create RNA-Seq Data Format (*.gct)
        tempTable = sortrows(R.expected_count,'geneName','ascend');
        tempTable = tempTable(:,[find(tempALocs),find(tempBLocs)]);
        
        % create gct file
        writetable(tempTable,sprintf('%s.gct',FileNameFull),...
        'Delimiter','\t','WriteVariableNames',1,'FileType','text')
        
        % add line 1 and 2 to .gct
        S = fileread(sprintf('%s.gct',FileNameFull));
        S = [sprintf('%i',height(tempTable)),char(9),...
            num2str(sum(tempALocs)+sum(tempBLocs)), char(10), S];
        S = ['#',FileName, char(10), S];
        FID = fopen(sprintf('%s.gct',FileNameFull), 'w');
        if FID == -1, error('Cannot open file %s', sprintf('%s.gct',FileNameFull)); end
        fwrite(FID, S, 'char');
        fclose(FID);
        
        % Create Phenotype Data Format (*.cls)
        fileID = fopen(sprintf('%s.cls',FileNameFull),'w');
        fprintf(fileID,'%i\t%i\t%i\n',num2str(sum(tempALocs)+sum(tempBLocs)),2,1);
        fprintf(fileID,'#%s\t%s\n',...
            sprintf('rnaseq_s%s_t%i',tempSampleA,tempTpA),...
            sprintf('rnaseq_s%s_t%i',tempSampleB,tempTpB));
        fprintf(fileID,sprintf('%s\n',num2str([repmat(1,1,sum(tempALocs)),repmat(2,1,sum(tempBLocs))])));
        fclose(fileID);
        
        % create names for combined table element
        nameAll{i} = sprintf('rnaseq_s%s_t%i',tempSampleA,tempTpA);
        nameAll{ii} = sprintf('rnaseq_s%s_t%i',tempSampleB,tempTpB);
    end
end

% add all diffExp comparisons to R structure
R.diffExpTable = cell2table(diffExpTable);

% name R.diffExpTable rows and columns
R.diffExpTable.Properties.VariableNames = nameAll;
R.diffExpTable.Properties.RowNames = nameAll;

end


function [R] = fdnDiffExpGsaa(dataInfo,R,gseaFlag,rnaseqPatternFlag)
%fdnDiffExpGsaa This function performs differential expression and
%creates files for GSAA
%
%   Input
%   dataInfo:   4DNvestigator structure that contains information on samples input
%   R:          4DNvestigator structure that contains RNA-seq data
%   gseaFlag:   flag for creating GSEA files
%   rnaseqPatternFlag:   flag for creating RNA-seq pattern heatmaps
%
%   Output 
%   R:          4DNvestigator structure that contains RNA-seq data with
%               differential expression table added
%
%   Version 2.0 (4/25/19)
%   Written by: Scott Ronquist
%   Contact:    scotronq@umich.edu
%   Created:    1/22/19
%   
%   Revision History:
%   v1.0 (1/22/19)
%   * fdnDiffExpGsaa.m created
%   v1.1 (03/15/19)
%   * formatted preamble
%   v2.0 (4/25/19)
%   * reorganized function
%       * new structure for diffExpTable (table in table)
%       * added gsea flag

%% Set up
if nargin < 3;gseaFlag = 0;end
if nargin < 4;rnaseqPatternFlag = 0;end

% determine number of samples
samplesAll = unique(dataInfo.sampleInfo.sample);

% create data table
diffExpCellSamp = cell(length(samplesAll),1);

%% Compare between all time points in a sample
for iS = 1:length(samplesAll)
    
    % determine time points within each sample
    sampleLocs = ismember(dataInfo.sampleInfo.sample,samplesAll{iS}) &...
        ismember(dataInfo.sampleInfo.dataType,'rnaseq');
    sampleTps = unique(dataInfo.sampleInfo.timePoint(sampleLocs));
    
    % create data structure for sample time point comparisons
    diffExpCellTp = cell(length(sampleTps));
    
    % loop through timepoints to compare between
    for iTp1 = 1:length(sampleTps)-1
        for iTp2 = iTp1+1:length(sampleTps)
            fprintf(['calculating Differential Expression for Sample "%s", ',...
                'TP %i vs TP %i...\n'],samplesAll{iS},sampleTps(iTp1),sampleTps(iTp2))
            
            % get sample time point locations
            sampleTpLoc1 = ismember(R.expected_count.Properties.VariableNames,...
                dataInfo.sampleInfo.name(ismember(dataInfo.sampleInfo.sample,samplesAll{iS}) &...
                ismember(dataInfo.sampleInfo.dataType,'rnaseq') &...
                ismember(dataInfo.sampleInfo.timePoint,sampleTps(iTp1))));
            sampleTpLoc2 = ismember(R.expected_count.Properties.VariableNames,...
                dataInfo.sampleInfo.name(ismember(dataInfo.sampleInfo.sample,samplesAll{iS}) &...
                ismember(dataInfo.sampleInfo.dataType,'rnaseq') &...
                ismember(dataInfo.sampleInfo.timePoint,sampleTps(iTp2))));
            
            % create geneTable
            tpmMeanTp1 = mean(R.TPM{:,sampleTpLoc1},2);
            tpmMeanTp2 = mean(R.TPM{:,sampleTpLoc2},2);
            meanBase = (tpmMeanTp1 + tpmMeanTp2) / 2;
            foldChange = tpmMeanTp1 ./ tpmMeanTp2;
            log2FC = log2(foldChange);
            
            geneTable = table(meanBase,tpmMeanTp1,tpmMeanTp2,foldChange,log2FC);
            geneTable = [R.TPM(:,1:6),geneTable];
            
            % DeSeq differential expression
            [pvalue,padj] = matlabNegbinDE(R.expected_count{:,sampleTpLoc1}, R.expected_count{:,sampleTpLoc2});
            geneTable.pvalue = pvalue;
            geneTable.padj = padj;
            
            % add to diffExpCellTp - level 2
            diffExpCellTp{iTp1,iTp2} = geneTable;
            
            %% GSAASeqSP analysis
            if gseaFlag
                % http://gsaa.unc.edu/userguide_gsaaseqsp.html
                % general info
                FileName = sprintf('%s_vs_%s',...
                    sprintf('rnaseq_s%s_t%i',samplesAll{iS},sampleTps(iTp1)),...
                    sprintf('rnaseq_s%s_t%i',samplesAll{iS},sampleTps(iTp2)));
                FileNameFull = fullfile(dataInfo.path.output,'data','gsaa',FileName);
                
                % Create RNA-Seq Data Format (*.gct)
                tempTable = sortrows(R.expected_count,'geneName','ascend');
                tempTable = tempTable(:,[find(sampleTpLoc1),find(sampleTpLoc2)]);
                
                % create gct file
                writetable(tempTable,sprintf('%s.gct',FileNameFull),...
                    'Delimiter','\t','WriteVariableNames',1,'FileType','text')
                
                % add line 1 and 2 to .gct
                S = fileread(sprintf('%s.gct',FileNameFull));
                S = [sprintf('%i',height(tempTable)),char(9),...
                    num2str(sum(sampleTpLoc1)+sum(sampleTpLoc2)), char(10), S];
                S = ['#',FileName, char(10), S];
                FID = fopen(sprintf('%s.gct',FileNameFull), 'w');
                if FID == -1, error('Cannot open file %s', sprintf('%s.gct',FileNameFull)); end
                fwrite(FID, S, 'char');
                fclose(FID);
                
                % Create Phenotype Data Format (*.cls)
                fileID = fopen(sprintf('%s.cls',FileNameFull),'w');
                fprintf(fileID,'%i\t%i\t%i\n',num2str(sum(sampleTpLoc1)+sum(sampleTpLoc2)),2,1);
                fprintf(fileID,'#%s\t%s\n',...
                    sprintf('rnaseq_s%s_t%i',samplesAll{iS},sampleTps(iTp1)),...
                    sprintf('rnaseq_s%s_t%i',samplesAll{iS},sampleTps(iTp2)));
                fprintf(fileID,sprintf('%s\n',num2str([repmat(1,1,sum(sampleTpLoc1)),repmat(2,1,sum(sampleTpLoc2))])));
                fclose(fileID);
            end
        end
    end
    
    % add to diffExpCellSamp - level 1
    diffExpCellSamp{iS,1} = cell2table(diffExpCellTp,...
        'RowNames',strcat('tp',cellstr(num2str(sampleTps))),...
        'VariableNames',strcat('tp',cellstr(num2str(sampleTps))));
    
    %% Determine gene expresssion patterns over time
    if rnaseqPatternFlag
        rnaseqPattern = diffExpCellSamp{iS,1}{1,2}{1}{:,8};
        for iTp = 2:length(sampleTps)
            rnaseqPattern = [rnaseqPattern diffExpCellSamp{iS,1}{1,iTp}{1}{:,9}];
        end
        
        % filter low expression genes and normalize over time
        rnaseqThresh = 1;
        rnaseqPatternNorm = rnaseqPattern;
        rnaseqPatternNorm(rnaseqPatternNorm<rnaseqThresh) = 0;
        rnaseqPatternNorm = normalize(rnaseqPatternNorm,2);
        
        % get expression patterns
        [C,IA,IC] = unique(diff(rnaseqPatternNorm,1,2)>0,'rows');
        
        % add NaN cluster
        IC(isnan(rnaseqPatternNorm(:,1))) = nan;
        
        % create red to green colormap
        n = 100;
        redGreen = colormap([linspace(1,0,n)', linspace(0,1,n)', zeros(n,1)] );  %// create colormap
        
        % sort expression patterns
        tempDataAll = [];
        tempDataNormAll = [];
        tempGenesAll = [];
        for iIC = 1:nanmax(IC)
            % get temp data
            tempData = rnaseqPattern(IC==iIC,:);
            tempDataNorm = rnaseqPatternNorm(IC==iIC,:);
            tempGenes = diffExpCellSamp{iS,1}{1,2}{1}{IC==iIC,4};
            
            % get total fold change
            tempFC = zeros(size(tempData,1),1);
            for iTp = 1:size(tempDataNorm,2)-1
                tempFC = tempFC + abs(log2(tempData(:,iTp)./tempData(:,iTp+1)));
            end
            
            % sort by total fold change
            [B,I] = sort(tempFC,'descend');
            tempData = tempData(I,:);
            tempDataNorm = tempDataNorm(I,:);
            tempGenes = tempGenes(I);
            
            % add to master table
            tempDataAll = [tempDataAll; tempData];
            tempDataNormAll = [tempDataNormAll; tempDataNorm];
            tempGenesAll = [tempGenesAll;tempGenes];
        end
        
        % plot image, all patterns together
        figure, imagesc(tempDataNormAll)
        colormap(redGreen), colorbar
        xticks(1:size(tempDataNorm,2)), xticklabels(sampleTps), xlabel('TP')
        ylabel('Normalized Gene Expression')
        title(sprintf('RNA-seq Patterns, TPM thresh=%i, Sample "%s"',rnaseqThresh,samplesAll{iS}))
    end
end

% add to main diffExpTable
R.diffExpTableTp = cell2table(diffExpCellSamp,...
    'RowNames',samplesAll,...
    'VariableNames',{'samples'});

%% Compare between all samples at a time points
% determine number of samples
tpAll = unique(dataInfo.sampleInfo.timePoint);

% create data table
diffExpCellTp= cell(length(tpAll),1);

%% Compare between all time points in a sample
for iTp = 1:length(tpAll)
    
    % determine time points within each sample
    tpLocs = ismember(dataInfo.sampleInfo.timePoint,tpAll(iTp)) &...
        ismember(dataInfo.sampleInfo.dataType,'rnaseq');
    tpSamples = unique(dataInfo.sampleInfo.sample(tpLocs));
    
    % create data structure for sample time point comparisons
    diffExpCellSample = cell(length(tpSamples));
    
    % loop through timepoints to compare between
    for iS1 = 1:length(tpSamples)-1
        for iS2 = iTp1+1:length(tpSamples)
            fprintf(['calculating Differential Expression for TP "%i", ',...
                'Sample "%s" vs Sample "%s"...\n'],tpAll(iTp),tpSamples{iS1},tpSamples{iS2})
            
            % get sample time point locations
            tpSamplLoc1 = ismember(R.expected_count.Properties.VariableNames,...
                dataInfo.sampleInfo.name(ismember(dataInfo.sampleInfo.timePoint,tpAll(iTp)) &...
                ismember(dataInfo.sampleInfo.dataType,'rnaseq') &...
                ismember(dataInfo.sampleInfo.sample,tpSamples{iS1})));
            tpSamplLoc2 = ismember(R.expected_count.Properties.VariableNames,...
                dataInfo.sampleInfo.name(ismember(dataInfo.sampleInfo.timePoint,tpAll(iTp)) &...
                ismember(dataInfo.sampleInfo.dataType,'rnaseq') &...
                ismember(dataInfo.sampleInfo.sample,tpSamples{iS2})));
            
            % create geneTable
            tpmMeanS1 = mean(R.TPM{:,tpSamplLoc1},2);
            tpmMeanS2 = mean(R.TPM{:,tpSamplLoc2},2);
            meanBase = (tpmMeanS1 + tpmMeanS2) / 2;
            foldChange = tpmMeanS1 ./ tpmMeanS2;
            log2FC = log2(foldChange);
            
            geneTable = table(meanBase,tpmMeanS1,tpmMeanS2,foldChange,log2FC);
            geneTable = [R.TPM(:,1:6),geneTable];
            
            % DeSeq differential expression
            [pvalue,padj] = matlabNegbinDE(R.expected_count{:,tpSamplLoc1}, R.expected_count{:,tpSamplLoc2});
            geneTable.pvalue = pvalue;
            geneTable.padj = padj;
            
            % add to diffExpCellTp - level 2
            diffExpCellTp{iS1,iS2} = geneTable;
            
            %% GSAASeqSP analysis
            if gseaFlag
                % http://gsaa.unc.edu/userguide_gsaaseqsp.html
                % general info
                FileName = sprintf('%s_vs_%s',...
                    sprintf('rnaseq_s%s_t%i',tpSamples{iS1},tpAll(iTp)),...
                    sprintf('rnaseq_s%s_t%i',tpSamples{iS2},tpAll(iTp)));
                FileNameFull = fullfile(dataInfo.path.output,'data','gsaa',FileName);
                
                % Create RNA-Seq Data Format (*.gct)
                tempTable = sortrows(R.expected_count,'geneName','ascend');
                tempTable = tempTable(:,[find(tpSamplLoc1),find(tpSamplLoc2)]);
                
                % create gct file
                writetable(tempTable,sprintf('%s.gct',FileNameFull),...
                    'Delimiter','\t','WriteVariableNames',1,'FileType','text')
                
                % add line 1 and 2 to .gct
                S = fileread(sprintf('%s.gct',FileNameFull));
                S = [sprintf('%i',height(tempTable)),char(9),...
                    num2str(sum(tpSamplLoc1)+sum(tpSamplLoc2)), char(10), S];
                S = ['#',FileName, char(10), S];
                FID = fopen(sprintf('%s.gct',FileNameFull), 'w');
                if FID == -1, error('Cannot open file %s', sprintf('%s.gct',FileNameFull)); end
                fwrite(FID, S, 'char');
                fclose(FID);
                
                % Create Phenotype Data Format (*.cls)
                fileID = fopen(sprintf('%s.cls',FileNameFull),'w');
                fprintf(fileID,'%i\t%i\t%i\n',num2str(sum(tpSamplLoc1)+sum(tpSamplLoc2)),2,1);
                fprintf(fileID,'#%s\t%s\n',...
                    sprintf('rnaseq_s%s_t%i',tpSamples{iS1},tpAll(iTp)),...
                    sprintf('rnaseq_s%s_t%i',tpSamples{iS2},tpAll(iTp)));
                fprintf(fileID,sprintf('%s\n',num2str([repmat(1,1,sum(tpSamplLoc1)),repmat(2,1,sum(tpSamplLoc2))])));
                fclose(fileID);
            end
        end
    end
    
    % add to diffExpCellSamp - level 1
    diffExpCellTp{iTp,1} = cell2table(diffExpCellSample,...
        'RowNames',tpSamples,...
        'VariableNames',tpSamples);
end

% add to main diffExpTable
R.diffExpTableSample = cell2table(diffExpCellTp,...
    'RowNames',strcat('tp',cellstr(num2str(tpAll))),...
    'VariableNames',{'tps'});

end


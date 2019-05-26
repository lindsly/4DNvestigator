classdef hic
    % hic is a class for storing Hi-C data
    %   hic stores raw data and normalization vectors stored in .hic files
    properties
        hicHeader = table();                    % table for hic header
        rawData  = {};                          % sparse matrices for raw data
        rawDataTable = cell2table(cell(0,2),...          
             'VariableNames',{'chr','res'});    % parameters run for raw data
        krVec = {};                             % vector for KR normalization
        oeVec = {};                             % vector for O/E normalization
        krEVec = {};                            % vector for KR O/E normalization
        badLocs = {};                           % vector for bad locs
    end
    
    methods
        function A = mat(obj,chrLoc1,chrLoc2,res,norm1d,norm3d)
            % This function extracts data from the objects stored within
            % this class
            
            %% Get default parameters
            % initialize base positioning
            basePos1 = [];
            basePos2 = [];
            
            % change string to char
            if isstring(chrLoc1); chrLoc1 = char(chrLoc1); end
            if isstring(chrLoc2); chrLoc2 = char(chrLoc2); end
            
            % Get base Chromosome
            if ischar(chrLoc1) && ischar(chrLoc2)
                if contains(chrLoc1,':') && contains(chrLoc1,':')
                    tempChrPos1 = strsplit(chrLoc1,':');
                    tempChrPos2 = strsplit(chrLoc2,':');
                    
                    baseChr1 = tempChrPos1{1};
                    baseChr2 = tempChrPos2{1};
                    
                    basePos1 = [str2num(tempChrPos1{2}) str2num(tempChrPos1{3})];
                    basePos2 = [str2num(tempChrPos2{2}) str2num(tempChrPos2{3})];
                elseif ~contains(chrLoc1,':') && ~contains(chrLoc1,':')
                    baseChr1 = chrLoc1;
                    baseChr2 = chrLoc2;
                else
                    error('chrLoc1 and chrLoc1 must input in the same format')
                end
            elseif ~ischar(chrLoc1) && ~ischar(chrLoc2)
                baseChr1 = num2str(chrLoc1);
                baseChr2 = num2str(chrLoc2);
            else
                error('chrLoc1 and chrLoc1 must input in the same format')
            end
            
            %% Error checks
            % Check if input resolution is valid
            if ~ismember(res,obj.hicHeader.BasePairdelimitedResolutions)
                error('user requested "res" is invalid')
            end
            
            % Check if input chr is valid
            if ~ismember(baseChr1,obj.hicHeader.Chromosomes.chr) ||...
                    ~ismember(baseChr2,obj.hicHeader.Chromosomes.chr)
                error('user requested "chr" is invalid')
            end
            
            %% Get index of rawDataTable
            if strcmp(baseChr1,baseChr2)
                rawDataTableIdx = strcmp(obj.rawDataTable.chr,baseChr1) &...
                    ismember(obj.rawDataTable.res,res);
            else 
                rawDataTableIdx = strcmp(obj.rawDataTable.chr,'ALL') &...
                    ismember(obj.rawDataTable.res,res);
            end
            
            if all(rawDataTableIdx)
                error('user requested "res" is not loaded for chr:%s and chr:%s',chrLoc1,chrLoc2)
            end
            
            %% Extract data
            % zeros pad, make square, make symmetric
            lengthA = length(obj.rawData{rawDataTableIdx});
            A = full(obj.rawData{rawDataTableIdx});
            A = padarray(A,lengthA-size(A),'post');
            A = A + triu(A,1)';
            
            % 1D norm
            switch lower(norm1d)
                case 'kr'
                    A = diag(obj.krVec{rawDataTableIdx}(1:lengthA).^(-1))*...
                        A*diag(obj.krVec{rawDataTableIdx}(1:lengthA).^(-1));
                case 'none'
                otherwise
                    error('user input "norm1d" is unrecognized')
            end
            
            % 3D norm
            switch lower(norm3d)
                case 'oe'
                    if lengthA > length(obj.krEVec{rawDataTableIdx})
                        A = A./toeplitz(padarray(obj.krEVec{rawDataTableIdx},...
                            lengthA-length(obj.krEVec{rawDataTableIdx}),'post'));
                    else
                        A = A./toeplitz(obj.krEVec{rawDataTableIdx}(1:lengthA));
                    end
                    
                case 'observed'
                otherwise
                    error('user input "norm3d" is unrecognized')
            end
            
            %% Extract sub-matrix
            if ~isempty(basePos1)
                A = A(ceil(basePos1(1)/res):ceil(basePos1(2)/res),...
                    ceil(basePos2(1)/res):ceil(basePos2(2)/res));
            end
        end
    end
end

function [stringLoc] = locInCellCell(cellCell,listOfStrings)
%locInCellCell find locations of a list in a cell of cell arrays
%   
%   Input
%   CellCell:       cell of cell arrays
%   listOfStrings: list of strings
%
%   Output
%   stringLoc:      Location of cell containing string from listOfStrings
%
%   Scott Ronquist, 1/22/19

%% Find strings in cell of cell array
stringLoc = cell(length(listOfStrings),1);
for i = 1:length(listOfStrings)
    stringLoc{i} = find(cellfun(@(x) any(strcmp(x,listOfStrings{i})),cellCell));
end

end


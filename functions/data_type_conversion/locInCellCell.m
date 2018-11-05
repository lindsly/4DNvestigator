function [stringLoc] = locInCellCell(CellCell,listOfStrings)
%locInCellCell find locations of a list in a cell of cell arrays
%   CellCell: cell of cell arrays
%   listOfInterest: list of strings
%
%   Scott Ronquist, 9/13/18

stringLoc = cell(length(listOfStrings),1);
for i = 1:length(listOfStrings)
    stringLoc{i} = find(cellfun(@(x) any(strcmp(x,listOfStrings{i})),CellCell));
end

end


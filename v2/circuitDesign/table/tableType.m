function string = tableType(table)
%TABLETYPE returns the type of a transistor (n or p)
%
%   string = tableType(TABLE) returns a string 'n' or 'p' depending on
%   whether the TABLE under consideration refers to an n-MOS or p-MOS
%   transistor. 
%
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%


if table.Info.Type == 'N'
  string = 'n';
else
  string = 'p';
end

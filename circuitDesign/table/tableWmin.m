function width = tableWmin(table)
%TABLEWMIN minimum channel width of a MOS transistor
%
%   WMIN = tableWmin(TABLE) returns the value of the minimum
%   channel width that is allowed for the transistor type that is specified by
%   the given TABLE. This value is usually smaller than the reference width.
%
%   EXAMPLE :
%
%     Wmin = tableWmin(N);
%
%   See also tableLmin, tableWcrit
%
%  (c) IMEC, 2004
%  IMEC confidential 
%


width = table.Info.Wmin;


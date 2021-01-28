function width = tableWcrit(table)
%TABLEWCRIT width of a MOS transistor below which narrow-channel effects play a role 
%
%   WCRIT = tableWcrit(TABLE) returns the value of the minimum
%   channel width that is allowed for the transistor type that is specified by
%   the given TABLE. This value is usually smaller than the reference width.
%
%   EXAMPLE :
%
%     wCrit = tableWcrit(N);
%
%   See also tableLmin, tableWref
%
%  (c) IMEC, 2004
%  IMEC confidential 
%


width = table.Info.Wcrit;


function length = tableLmin(table)
%TABLELMIN minimum channel length of a MOS transistor
%
%   LMIN = tableLmin(TABLE) returns the value of the minimum
%   channel length that is allowed for the transistor type that is specified by
%   the given TABLE of intrinsic operating point parameters of a MOS transistor.
%
%   EXAMPLE :
%
%     lmin = tableLmin(N);
%
%   See also tableWmin, tableWcrit
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

length = table.Info.Lmin;


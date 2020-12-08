function w = tableWref(table)
%TABLEWREF gets the MOS transistor width with which the table of that MOS 
%is computed.
%
%   W = tableWref(TABLE) returns the width in micrometer of the transistors with
%   which the values in the given TABLE have been computed.
%
%   EXAMPLE :
%
%      width = tableWref(N);
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

w = table.Info.Wref;




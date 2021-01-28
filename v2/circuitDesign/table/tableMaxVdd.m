function value = tableMaxVdd(table)
%TABLEMAXVDD retrieving the maximally allowed VDD for a MOS transistor
%
%   VALUE = tableMaxVdd(TABLE) returns the value of the maximally allowed VDD for
%   a MOS transistor that is specified by a given TABLE of intrinsic
%   operating point parameters of a MOS transistor.
%
%   EXAMPLE :
%
%     vdd = tableMaxVdd(N);
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

value = table.Info.vddMax;


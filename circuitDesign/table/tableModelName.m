function string = tableModelName(table)
%TABLEMODELNAME model name of MOS transistor with which the table values 
% are computed.
%
%   STRING = tableModelName(TABLE) returns the name of the MOS model with which
%   the values in the given TABLE have been computed.
%
%   EXAMPLE :
%
%     modelname = tableModelName(N);
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

string = table.Info.modelName;

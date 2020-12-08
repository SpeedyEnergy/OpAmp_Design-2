function list = tableWInverseParamNames(table)
%TABLEWINVERSEPARAMNAMES returns operating point parameters of a MOS that are 
% assumed to be inversely proportional to the transistor width. Currently, this is
% only the power spectral density of the 1/f noise current source between
% source and drain ("di2_fn").
%
%   stringArray = tableWInverseParamNames(TABLE) returns a cell array of
%   strings with the names of operating point parameters stored in the given
%   TABLE, and that - above a critical width - are proportional to 
%   (transistor width)^(-1)
%
%   EXAMPLE :
%
%    stringArray = tableWInverseParamNames(N);
%
%   See also tableWcrit, tableWIndepParamNames
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

list = table.Info.WInverseParams;

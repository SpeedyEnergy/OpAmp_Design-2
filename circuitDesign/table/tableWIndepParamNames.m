function list = tableWIndepParamNames(table)
%TABLEWINDEPPARAMNAMES returns operating point parameters of a MOS that are 
%assumed to be width independent. 
%
%   stringArray = tableWIndepParamNames(TABLE) returns a cell array of
%   strings with the names of operating point parameters stored in the given
%   TABLE, and that are independent of transistor width above a critical width.
%
%   EXAMPLE :
%
%    stringArray = tableWIndepParamNames(N);
%
%   See also tableWcrit
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

list = table.Info.WIndepParams;

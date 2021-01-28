function yesOrNo = isModelParam(string, table)
%ISMODELPARAM returns true if a given string is a model parameter of a ...
%transistor characterized by a table.
%
%   yesOrNo = isModelParam(STRING, TABLE) returns 1 if the given STRING is
%   defined as a model parameter for the type of transistor defined by the
%   given TABLE 
%
%   EXAMPLE :
%
%    yOrN = isModelParam('pbsw', N)
%
%

yesOrNo = ismember(string, tableModelParamNames(table));

function result = tableModelParamNames(table)
%TABLEMODELPARAMNAMES list of all names of model parameters 
%
%   RESULT = tableModelParamNames(TABLE) returns a cell array of strings. This cell
%   array contains the names of all model parameters that have been defined in
%   the given TABLE of a MOS transistor.
%
%   EXAMPLE :
%
%    modelParamNames = tableModelParamNames(N);
%
%   See also tableModelParam, tableModelName
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%


result = table.Info.modelParams;


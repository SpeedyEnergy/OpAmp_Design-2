function value = tableModelParam(parameter, table)
%TABLEMODELPARAM getting the value of a model parameter from a MOS table
%
%   VALUE = tableModelParam(PARAMETER, TABLE) retrieves from the data of a given
%   MOS table TABLE the numerical value of a model parameter that is specified
%   with the string PARAMETER. No unit conversion is performed, which means
%   that modelparameters that are not specified in S.I. units are returned in
%   non-S.I. units.
%
%   EXAMPLE :
%
%      M1.cjsw = tableModelParam('cjsw', N);
%
%    See also tableModelParamNames, tableModelName
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%


if not(ismember(parameter, tableModelParamNames(table)))
  error('the specified model parameter %s is not specified in the table\n', parameter);
end

value = table.Model.(parameter);


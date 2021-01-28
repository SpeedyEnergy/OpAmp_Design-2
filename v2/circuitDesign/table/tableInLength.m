function result = tableInLength(inputVar, table)
%TABLEINLENGTH getting the number of values of an input variable of a MOS table.
%
%    N = tableInLength(INPUTVAR, TABLE) returns the number
%    of values of the input variable INPUTVAR for a given TABLE of a MOS
%    transistor.
%
%    EXAMPLE :
%
%      nVgs = tableInLength('vgs', N);
%
%    See also tableInInit, tableInFinal, tableInArray, tableInValue, tableInStep,
%    tableInLength.
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

if not(ismember(inputVar, tableInNames(table)))
  error('the specified input argument is not an input variable');
end

result = table.Input.(inputVar).nInputs;


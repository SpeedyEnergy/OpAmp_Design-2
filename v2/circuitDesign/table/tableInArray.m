function result = tableInArray(inputVar, table)
%TABLEINARRAY getting the array of values of an input variable of a MOS table.
%
%    VALUE = tableInArray(INPUTVAR, TABLE) returns the array
%    of values of the input variable INPUTVAR for a given TABLE of a MOS
%    transistor.
%
%    EXAMPLE :
%
%      vgsArray = tableInArray('vgs', N);
%
%    See also tableInFinal, tableInInit, tableInValue, tableInStep,
%    tableInLength.
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

if not(ismember(inputVar, tableInNames(table)))
  error('the specified input argument is not an input variable');
end

result = table.Input.(inputVar).array;

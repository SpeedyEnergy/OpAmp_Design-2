function result = tableInInit(inputVar, table)
%TABLEININIT getting the initial value of an input variable of a MOS table.
%
%    VALUE = tableInInit(INPUTVAR, TABLE) returns the first element of the array
%    of values of the input variable INPUTVAR for a given TABLE of a MOS
%    transistor.
%
%    EXAMPLE :
%
%      vgsInit = tableInInit('vgs', N);
%
%    See also tableInFinal, tableInArray, tableInValue, tableInStep,
%    tableInLength.
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

if not(ismember(inputVar, tableInNames(table)))
  error('the specified input argument is not an input variable');
end

result = table.Input.(inputVar).init;

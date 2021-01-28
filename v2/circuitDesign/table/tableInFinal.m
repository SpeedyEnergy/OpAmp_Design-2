function result = tableInFinal(inputVar, table)
%TABLEINFINAL getting the final value of an input variable of a MOS table.
%
%    VALUE = tableInFinal(INPUTVAR, TABLE) returns the last element of the array
%    of values of the input variable INPUTVAR for a given TABLE of a MOS
%    transistor.
%
%    EXAMPLE :
%
%      vgsFinal = tableInFinal('vgs', N);
%
%    See also tableInInit, tableInArray, tableInValue, tableInStep,
%    tableInLength.
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

if not(ismember(inputVar, tableInNames(table)))
  error('the specified input argument is not an input variable');
end

result = table.Input.(inputVar).final;

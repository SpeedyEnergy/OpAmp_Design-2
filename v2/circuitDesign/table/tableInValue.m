function result = tableInValue(inputVar, index, table)
%TABLEINVALUE gets a particular value from the array of values 
%of an input variable of a MOS table of intrinsic operating point parameters.
%
%    VALUE = tableInValue(INPUTVAR, INDEX, TABLE) returns the element with
%    index INDEX in the array of values of the input variable INPUTVAR in a
%    given TABLE of a MOS transistor.
%
%    EXAMPLE :
%
%      vgsInter = tableInValue('vgs', 5, N);
%
%    See also tableInFinal, tableInInit, tableInArray, tableInStep,
%    tableInLength.
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

if not(ismember(inputVar, tableInNames(table)))
  error('the specified input argument is not an input variable');
end

if not(isnumeric(index)) | (index < 1) | (index > tableInLength(inputVar, ...
      table))
  error('The specified index is not in the correct range');
end

result = table.Input.(inputVar).array(index);

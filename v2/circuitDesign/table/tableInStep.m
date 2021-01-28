function result = tableInStep(inputVar, table)
%TABLEINSTEP step between two adjacent values of an input variable of a MOS table.
%
%    VALUE = tableInStep(INPUTVAR, TABLE) returns the value of the step
%    of the input variable INPUTVAR for a given TABLE of a MOS
%    transistor. This step can be positive or negative. When the step is not
%    constant, then the function returns NaN.
%
%    EXAMPLES :
%
%      vgsStep = tableInStep('vgs', N);
%      lgStep = tableInStep('lg', N);
%
%       This second example will return NaN.
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

result = table.Input.(inputVar).step;

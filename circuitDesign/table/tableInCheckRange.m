function tableInCheckRange(inputVar, value, table);
%TABLEINCHECKRANGE check whether an input parameter is inside the validity range
%
%   tableInCheckRange(INPUTVAR, VALUE, TABLE) tests whether the VALUE of a given
%   input parameter INPUTVAR (specified as a string) is inside the range of
%   validity, i.e. the range that is used to construct the given TABLE for a
%   MOS transistor.
%   The function signals an error if the VALUE is not inside the range of
%   validity.
%
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

minVal = min(tableInInit(inputVar, table), tableInFinal(inputVar, table));
maxVal = max(tableInInit(inputVar, table), tableInFinal(inputVar, table));
if (value < minVal - eps) | (value > maxVal + eps)
    error('value %g is out of range of %s values', value, inputVar);
end


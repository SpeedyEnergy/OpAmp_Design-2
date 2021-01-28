function list = cirDefaultValues(circuit, eltype)
%CIRDEFAULTVALUES lists the default values of the fields that are initialized 
%when a circuit element is defined.
%
%    LIST = cirDefaultValues(CIRCUIT, ELTYPE) returns a cell array with the
%    default values to which the fields of a circuit element of type ELTYPE
%    are initialized when this circuit element is defined with cirElementDef.
%
%    See also cirInit, cirElementDef, cirDefaultNames. 
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

list = circuit.defaults.(eltype).fieldValues;

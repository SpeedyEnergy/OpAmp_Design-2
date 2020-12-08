function list = cirDefaultNames(circuit, eltype)
%CIRDEFAULTNAMES lists the fields names that are initialized when a circuit 
%element is defined.
%
%    LIST = cirDefaultNames(CIRCUIT, ELTYPE) returns a cell array of
%    strings that are the fieldnames of a circuit element of type ELTYPE
%    when this circuit element is defined with cirElementDef.
%
%    See also cirInit, cirElementDef, cirDefaultValues. 
%
%  (c) IMEC, 2004
%  IMEC confidential 
%


list = circuit.defaults.(eltype).fieldNames;


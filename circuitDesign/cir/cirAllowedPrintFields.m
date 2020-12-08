function list = cirAllowedPrintFields(circuit, eltype)
%CIRALLOWEDPRINTFIELDS returns a list names of the fields of a circuit
%    element of type ELTYPE, that can be dumped into/modified in a circuit
%    simulation file. 
%
%    See also cirInit
%
%    EXAMPLE:
%
%    cirAllowedPrintFields(circuit, 'nmos')
%
%    This returns the list {'w', 'lg', 'ad', 'as', 'pd', 'ps', ...
%    'nFingers', 'lsos', 'lsod', 'lsogs', 'lsogd', 'mult'}
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

cirCheckEltype(eltype);
list = circuit.allowedPrintFields.(eltype);

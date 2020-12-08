function elementList = cirElementOfType(circuit, eltype)
%CIRELEMENTOFTYPE returns a list of elements of a given type.
%
%    ELEMENTLIST = cirElementOfType(CIRCUIT, ELTYPE) returns a cell array of
%    circuit elements of a given type ELTYPE (specified as a string), that
%    have been defined in the given CIRCUIT.
%
%    See also cirInit, cirDefaultNames, cirDefaultValues.
%    
%    EXAMPLES:
%
%    nmosList = cirElementOfType(cir, 'nmos')
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%


elementListLength = 0;
fieldList = fieldnames(circuit);
for i = 1:length(fieldList)
    if (isfield(circuit.(fieldList{i}), 'eltype') & strcmp(eltype, circuit.(fieldList{i}).eltype))
        elementListLength = elementListLength + 1;
        elementList{elementListLength} = circuit.(fieldList{i});
    end
end

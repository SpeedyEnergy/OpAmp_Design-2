function destination = cirElementCopy(source, destination)
%CIRELEMENTCOPY copies a circuit element to another one.
%
%    DESTINATION = cirElementCopy(SOURCE, DESTINATION) copies a circuit
%    element SOURCE to another circuit element DESTINATION. By this copy
%    operation, all fields except for the name field are copied from the
%    SOURCE element to the DESTINATION element. It is required that the
%    element types of SOURCE and DESTINATION correspond. 
%
%
%    EXAMPLES:
%
%    Mn1b = cirElementCopy(Mn1a, Mn1b);
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%


if not(strcmp(source.eltype, destination.eltype))
    error('Copying circuit elements of different types is not allowed');
end

destinationName = destination.name;
destination = source;
destination.name = destinationName;


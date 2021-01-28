function filename = cirSimulSkelFile(circuit)
%CIRSIMULSKELFILE retrieves the circuit simulation file
%
%   FILENAME = cirSimulSkelFile(CIRCUIT) returns a string that contains for the
%   given CIRCUIT the name of the file that contains the part of the circuit
%   simulation input file that is fixed (i.e. contains no elements that are
%   determined with Matlab).
%
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

filename = circuit.simulSkeletonFile;


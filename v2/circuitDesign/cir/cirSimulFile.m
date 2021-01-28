function filename = cirSimulFile(circuit)
%CIRSIMULFILE retrieves the circuit simulation file
%
%   FILENAME = cirSimulFile(CIRCUIT) returns a string that contains for the
%   given CIRCUIT the name of the file that will be used to verify the design
%   with Matlab of this CIRCUIT.
%
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

filename = circuit.simulFile;


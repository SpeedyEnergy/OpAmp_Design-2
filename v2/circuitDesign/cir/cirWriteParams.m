function cirWriteParams(circuit, simulator)
%CIRWRITEPARAMS writes parameters that have been determined in Matlab to a
%circuit simulation file for a given SIMULATOR.
% 
%   Supported simulators are smartspice, hspice, spectre
%
%
%   For smartspice we put everything into one file in order to avoid .include
%   statements, which seem to be buggy in smartspice. Therefore, we dump the
%   parameters in a temporary file, and in the end we combine this file and
%   the skeleton file into the simulation file.
%
%   For HSPICE, the simulation file is just a set of lines with .param
%   statements. The skeleton file is not used here.
%
%   For spectre, we replace in a circuit netlist the values of the circuit
%   elements in the top circuit and in the subcircuits. The skeleton file is
%   not used neither here.
%
%  See also cirWriteParamsSpice and cirWriteParamsSpectre
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

switch simulator
  case {'smartspice', 'hspice'}
    cirWriteParamsSpice(circuit, simulator);
  case 'spectre'
    cirWriteParamsSpectre(circuit, simulator);
  otherwise
    error('specified simulator not supported');
end

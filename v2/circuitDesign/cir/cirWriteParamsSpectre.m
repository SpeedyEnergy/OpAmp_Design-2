function parStruct = cirWriteParamsSpectre(varargin)
%CIRWRITEPARAMS writes parameters that have been determined in Matlab to a
%   circuit simulation file for spectre. We replace
%   in a circuit netlist the values of the circuit elements in the top circuit
%   and in the subcircuits. We use the function scsSetPar which substitutes
%   in the simulation file circuit.simulFile, using a structure of parameters. 
%   This structure is constructed in this function cirWriteParamsSpectre. 
%   The structure can be saved in a .mat file, which can be given as a second argument.
%   This function also strips unnecessary lines from a spectre netlist,
%   that can be generated automatically by the Cadence Analog Design
%   Environment. For the removal of unnecessary lines, scsStripNetlist2 is
%   used, not scsStripNetlist.
%
%   See also scsSetPar, scsStripNetlist2
%
%  EXAMPLE:
% 
%  cirWriteParamsSpectre(circuit)
%  Here the substitution is made without saving the parameter structure.
%  The substitution is done in the file circuit.simulFile.
%
%  cirWriteParamsSpectre(circuit, 'parStruct.mat')
%  Here the parameter structure is saved to a file parStruct.mat
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

debug = 0;

switch nargin
  case 1
    saveParStruct = 0;
    circuit = varargin{1};
  case 2
    circuit = varargin{1};
    parStructFileName = varargin{2};
  otherwise
    error('wrong number of arguments');
end

if not(isfield(circuit, 'simulFile'))
    error('The circuit does not have a simulation file (field "simulFile")');
end
scsStripNetlist2(circuit.simulFile);
parStruct.dummy = 1;
eltypes = cirSupportedEltypes;
for i = 1:length(eltypes)
  nElementsField = strcat('n', eltypes{i});
  % this is the name of the field of the given circuit that contains the number
  % of circuit elements of the given eltype
  elementsListField = strcat(eltypes{i}, 'List'); 
  % this is the name of the field of the given circuit that contains the list
  % of all circuit elements of the given eltype
  
  for j = 1:circuit.(nElementsField)
    element = circuit.(circuit.(elementsListField){j});
    if debug
      fprintf(1, 'calling addSpectreParam with element %s\n', element.name);
    end
    parStruct = addSpectreParam(element, parStruct); 
  end
end
parStruct = rmfield(parStruct, 'dummy')


fprintf(1, 'saving the substituted parameters in file %s\n', ...
    parStructFileName); 
save parStructFileName parStruct; 
origPar = scsSetPar(circuit.simulFile, parStruct);






function parStruct = addSpectreParam(element, parStruct)
debug = 0;

switch element.eltype
  case 'subckt'
    parStruct.(element.name).dummy = 1;
    eltypes = cirSupportedEltypes;
    for i = 1:length(eltypes)
      nElementsField = strcat('n', eltypes{i});
      % this is the name of the field of the given circuit that contains
      % the number of circuit elements of the given eltype
      elementsListField = strcat(eltypes{i}, 'List'); 
      % this is the name of the field of the given circuit that contains
      % the list of all circuit elements of the given eltype
      for j = 1:element.(nElementsField)
	if debug
	  fprintf(1, 'element.(elementsListField){j} is \n\n');
	  element.(elementsListField){j}
	  fprintf(1, '\n\n\n');
	end
	parStruct.(element.name) = ...
	    addSpectreParam(element.(element.(elementsListField){j}), ...
	    parStruct.(element.name));
      end
    end
    parStruct.(element.name) = rmfield(parStruct.(element.name), 'dummy');
  otherwise
    for j = 1:length(element.paramFields)
      parStruct.(element.name).(element.paramPrintNames{j}) = ...
	  element.(element.paramFields{j});  
    end
end

    

	      
      

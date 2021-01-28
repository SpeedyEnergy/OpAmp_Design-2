function circuit = cirElementDef(varargin)
%CIRELEMENTDEF definition of a circuit element.
%
%    CIRCUIT = cirElementDef(CIRCUIT, ELTYPE, ELEMENTNAME) adds a circuit
%    element to the datastructure of the given CIRCUIT. This function changes
%    the datastructure of the given circuit, which is returned. The circuit
%    element needs to be specified by a name ELEMENTNAME, which should be a
%    string and by its element type ELTYPE, which should also be a string.
%    The circuit element is added as a field to the circuit, such as
%    CIRCUIT.ELEMENTNAME. This is in turn a structure, which will have
%    different fields depending on the element type. Two fields are already
%    defined in this function, namely CIRCUIT.ELEMENTNAME.eltype and
%    CIRCUIT.ELEMENTNAME.name, which are filled with the strings ELTYPE and
%    ELEMENTNAME, respectively.
%    Further, the field nElements of the given CIRCUIT is increased with one,
%    the field CIRCUIT.(ELTYPELIST) is augmented with the ELEMENTNAME, and the
%    field CIRCUIT.n(ELTYPE) is augmented with one. For example, for an n-MOS
%    transistor Mn1, eltype is 'nmos', and the field CIRCUIT.nmosList is
%    extended with the name 'Mn1', and the field CIRCUIT.nnmos is augmented
%    with one. 
%
%    For the circuit element itself, the field "parent" is defined and it is
%    equal to the name of the circuit. 
%
%    In addition to the fields eltype and name, some other fields, that are 
%    specific for each element type (see cirDefaultNames(CIRCUIT, ELTYPE)), 
%    are already defined and initialized with default values.
%    These default values can be overridden during initialization
%    by specifying extra arguments to cirElementDef in addition to the three
%    arguments CIRCUIT, ELTYPE and ELEMENTNAME. If no extra arguments are
%    specified, then default values are used to initialize the fields of 
%    cirDefaultNames(CIRCUIT, ELTYPE). Else, when default values need to be 
%    overridden, the syntax is as follows:
%    CIRCUIT = cirElementDef(CIRCUIT, ELTYPE, ELEMENTNAME, 'fieldname', value,
%    'fieldname', value, ...)
%    This means that the fieldname (string) has to be specified as an argument,
%    immediately followed by its value as the next argument. 
%
%    See also cirElementsCheckIn, cirDefaultNames, cirDefaultValues.
%    
%    EXAMPLES:
%
%    cir = cirElementDef(cir, 'nmos', 'Mn1')
%
%         This example defines n-MOS transistor Mn1 and takes default values
%         for the fields lg, w, table, paramFields, paramPrintNames, lsos,
%         lsod, mult, lsogs, lsogd, mult.
%
%    cir = cirElementDef(cir, 'nmos', 'Mn1', 'lg', 2e-6)
%
%        This example defines n-MOS transistor Mn1. The channel length lg is
%        initialized to 2 micrometer, while the other fields w, table, 
%        paramFields, paramPrintNames, lsos, lsod, mult, lsogs, lsogd, mult 
%        are initialized to their default value.
%
%   cir = cirElementDef(cir, 'nmos', 'Mn3', 'paramFields', {nFingers})
%
%       This example defines n-MOS transistor Mn3. The field
%       "paramField" is initialized with {nFingers}. This means that only
%       the value field "nFingers" will be printed to the circuit simulation
%       file. The other fields are initialized with their default value.
%
%  cir = cirElementDef(cir, 'nmos', Mn4, 'paramFields', {gm, lg, w})
%       This will give an error, since "gm" is not a parameter that can be
%       specified as one that will be printed to the circuit simulation file.
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%




if nargin < 3
  error('Too few arguments');
end

if not(mod(nargin,2))
  error('Incorrect number of arguments');
end


circuit = varargin{1};
eltype = char(varargin{2});
elementName = char(varargin{3});

if isfield(circuit, elementName)
  error('The specified name already exists');
end

if not(ischar(eltype))
  error('The specified circuit element type is not a string');
end

circuit.nElements = circuit.nElements + 1;

eltypeListField = strcat(eltype, 'List'); 
% this is the name of the field of the given circuit that contains the list of
% all circuit elements of the given eltype

eltypeNumberField = strcat('n', eltype);
% this is the name of the field of the given circuit that contains the number
% of circuit elements of the given eltype

circuit.(elementName).eltype = eltype;
circuit.(elementName).name = elementName;
circuit.(elementName).parent = circuit.name;
circuit.(eltypeNumberField) = circuit.(eltypeNumberField) + 1;
circuit.(eltypeListField){circuit.(eltypeNumberField)} = elementName;

% rest of the fields is first filled with default values (not for subcircuits):
if not(strcmp(eltype, 'subckt')) 
  elementFields = cirDefaultNames(circuit, eltype);
  elementFieldValues = cirDefaultValues(circuit, eltype);
  for i = 1:length(elementFields)
    circuit.(elementName).(char(elementFields{i})) = elementFieldValues{i};
  end
  
  % next the defaults are overridden when the field under consideration and its
  % new value have been specified as an argument to this function:
  for i=4:2:length(varargin)
    %fprintf(1, 'varargin{%d} is \n', i);
    %varargin{i}
    %fprintf(1, '\n\n');
    
    if not(ischar(varargin{i}))
      error('Argument is not a string');
    end
    if not(ismember(varargin{i}, elementFields))
      error('Argument unknown for initialization of this circuit element');
    end
    circuit.(elementName).(varargin{i}) = varargin{i+1};
  end
end


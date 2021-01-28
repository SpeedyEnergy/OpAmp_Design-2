function circuit = cirInit(name, circuitTitle, hierarchy, elementList, ...
    varargin)  
%CIRINIT initialization of a circuit
%
% circuit = cirInit(NAME, TITLE, HIERARCHY, ELEMENTLIST, ...
%   SPECS (optional), DESIGNCHOICES (optional), DESIGNKITNAME (optional),
%   NMOSTABLE (optional), PMOSTABLE (optional), SIMULATOR (optional), ...
%   SIMULFILE (optional), PARAMSIMULSKELFILE (optional)) 
%   initializes a datastructure for a circuit. This new datastructure is
%   returned. If no hierarchy is built in the circuit, then the arument
%   HIERARCHY must have the value 'top'. In a hierarchical circuit, this
%   argument must have this value as well. If a subcircuit is initialized, then
%   the argument HIERARCHY must have the value 'sub'.
%
%   The datastructure of a circuit is a structure with the following fields:
%
%     1. name: a string, which is specified as the argument NAME of this
%        initialization function.
%     
%     2. title: a title string, which is specified as the argument TITLE 
%        of this initialization function.
%
%     3. simulFile: a string that represents the name of the file to which the
%        values will be written for the different circuit parameters that are
%        determined with Matlab. This string is specified as the argument
%        SIMULFILE to this initialization function.
%
%     4. simulSkeletonFile: this field is only used for circuit simulations in 
%        SmartSpice. It is a string representing the name of the file that
%        contains the part of a circuit simulation file that does not contain
%        any parameters that are determined in Matlab.  Therefore, it is an optional
%        argument SIMULSKELFILE, that does not need to be specified when other 
%        simulators are used. 
%     
%     5. defaults: a structure that contains default values for the different
%        circuit elements that are supported. These default values are used when a
%        circuit element is initialized. For a nMOS and pMOS transistor,
%        default values are used for 
%        - channel length (field "lg")
%        - channel width (field "w")
%        - table of intrinsic operating point parameters (field "table")
%        - length of drain/source diffusion areas not between two
%         poly stripes (fields "lsod" and "lsos", respectively)
%        - length of drain/source diffusion areas between two poly stripes
%         (fields "lsogd" and "lsogs", respectively)
%        - the amount of fingers (field "nFingers")
%        - the multiplicity (field "mult")
%        - the parameters that will be written
%          to a circuit simulation file (field "paramFields"). This is a cell
%          array of strings, each of them being a field of a MOS transistor
%          structure (e.g. "lg"). The value of these fields will be determined
%          in the Matlab design plan.
%        - the printnames that will be used for the parameters in the field
%          paramFields of the transistor (field "paramPrintNames"). This is also
%          a cell array of strings. This cell array has the same length as the
%          cell array "paramFields". The printnames are used for the names of
%          the parameters that are written into the circuit simulation
%          file. E.g. the printname of the length of a MOS transistor could be
%          'l'. In that case, a parameter would be made of the form
%          '<mos><glue string>l', e.g. m1_l
%     For a voltage source, current source, capacitor, inductor and resistor, 
%     default values are used for 
%       - value (field "val")
%       - field "paramFields" (similar meaning as for a MOS transistor)
%       - field "paramPrintNames" (similar meaning as for a MOS transistor)
%
%     6. allowedPrintFields: an array of strings (cell array) that contains the
%     names of the fields that can be passed to/changed in a circuit simulation
%     file. For example, for a MOS transistor this list is
%     {'w', 'lg', 'ad', 'as', 'pd', 'ps', 'nFingers', 'lsos', 'lsod', ...
%      'lsogs', 'lsogd', 'mult', 'nrs', 'nrd'}.
%     For a resistor, capacitor, inductor, voltage source and current source,
%     this list contains one element, namely 'val'.
%
%     7. nElements: number of circuit elements that already have been defined.
%
%     8. fields that each represent a circuit element that belongs to the
%        circuit. At initialization there is no such field. During the Matlab
%        program these fields are filled. E.g., when an n-MOS transistor m1
%        is defined (see function cirElementDef) then there will be a field
%        "m1" in the circuit.
%        The following types of elements are supported: nmos, pmos, resistor,
%        capacitor, voltage source, current source (see function
%        cirCheckEltype). Note that not every circuit element that appears in
%        the netlist of the circuit, needs to be defined in Matlab. Only those
%        elements need to be defined for which some parameters are computed in
%        Matlab. For example, if the power supply voltage is fixed, then the
%        corresponding voltage source does not need to be added to the
%        datastructure. 
%    
%     9. fields that represent for each circuit element type a list of the
%        elements of that type and the length of that list:
%        - nmosList and nnmos
%        - pmosList and npmos
%        - resList and nres, capList and ncap, indList and nind, vList and
%          nv, iList and ni, subcktList and nsubckt
%
%     10. designkit
%
%     11. if the argument SPECS has been specified as an argument that differs
%         from 0 or from an empty string, then the circuit that is initialized
%         gets a field for its specifications. This field is meant to be a
%         structure. This is useful to store with the circuit sizing.
%         Example: 
%              if a structure "spec" has the following fields
%                  spec.gain = 5
%                  spec.currentConsumption = 1e-3;
%              then the inclusion of "spec" in the argument list of cirInit to
%              initialize a circuit with name "amp", will result in the field 
%              amp.spec and the value will be
%                  amp.spec.gain = 5
%                  and
%                  amp.spec.currentConsumption = 1e-3;
%
%    12. if the argument DESIGNCHOICES has been specified as an argument that differs
%         from 0 or from an empty string, then the circuit that is initialized
%         gets a field for its design choices. This field is meant to be a
%         structure. This is useful to store with the circuit sizing.
%         Example: 
%              if a structure "choice" has the following fields
%                  choice.Mn1.lg = 1e-6
%                  choice.Ibias.val = 1e-6
%              then the inclusion of "choice" in the argument list of cirInit to
%              initialize a circuit with name "amp", will result in the field 
%              amp.choice and the value will be
%                  amp.choice.Mn1.lg = 1e-6
%                  and
%                  amp.choice.Ibias.val = 1e-6
%             
%
% See also cirSubcktInit, cirElementDef, cirElementOfType, cirCheckElType,
% cirSupportedEltypes, cirDesignkitInit 
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

debug = 0;
if debug
  fprintf(1, 'argument 2 = %s\n', inputname(2))
  fprintf(1, 'variable argument 1 = \n');
  disp(varargin{1});
  fprintf(1, 'inputname(4) is %s\n', inputname(4));
  fprintf(1, 'elementList is\n');
  disp(elementList);
  fprintf(1, '\n');
  if nargin > 4
    fprintf(1, 'inputname(5) is %s\n', inputname(5));
  end
  if nargin > 5
    fprintf(1, 'inputname(6) is %s\n', inputname(6));
    fprintf(1, 'variable argument 2 = \n');
    disp(varargin{2});
  end
  if nargin > 6
    fprintf(1, 'inputname(7) is %s\n', inputname(7));
    fprintf(1, 'variable argument 3 = \n');
    disp(varargin{3});
  end
  if nargin > 7
    fprintf(1, 'inputname(8) is %s\n', inputname(8));
    fprintf(1, 'variable argument 4 = \n');
    disp(varargin{4});
  end
  if nargin > 8
    fprintf(1, 'inputname(9) is %s\n', inputname(9));
  end
  if nargin > 9
    fprintf(1, 'inputname(10) is %s\n', inputname(10));
  end
  if nargin > 10
    fprintf(1, 'inputname(11) is %s\n', inputname(11));
  end
  if nargin > 11
    fprintf(1, 'inputname(12) is %s\n', inputname(12));
  end
  
  fprintf(1, '\n\n');
end


if (nargin < 4) | (nargin > 12)
  error('Wrong number of arguments');
end

circuit.name = name;
circuit.title = circuitTitle;
if nargin > 4
  specs = varargin{1};
  if isstruct(specs)
    circuit.(inputname(5)) = specs;
  end
end

if nargin > 5
  choices = varargin{2};
  if isstruct(choices)
    circuit.(inputname(6)) = choices;
  end
end

if nargin > 6
  designkitName = varargin{3};
end

%%%%%%%%%%%%%%%
if nargin > 7
  nmosTable = varargin{4};
  nmosTableName = inputname(8);
end
if nargin > 8
  pmosTable = varargin{5};
  pmosTableName = inputname(9);
end
if nargin > 9
  simulator = varargin{6};
end
if nargin > 10
  circuit.simulFile = varargin{7};
end
if nargin == 12
  circuit.simulSkeletonFile = varargin{8};
end


circuit.allowedPrintFields.nmos = {'w', 'lg', 'ad', 'as', 'pd', 'ps', ...
    'nFingers', 'lsos', 'lsod', 'lsogs', 'lsogd', 'mult', 'nrs', ...
    'nrd', 'nFins'};
circuit.allowedPrintFields.pmos = {'w', 'lg', 'ad', 'as', 'pd', 'ps', ...
    'nFingers', 'lsos', 'lsod', 'lsogs', 'lsogd', 'mult', 'nrs', ...
    'nrd', 'nFins'};

circuit.allowedPrintFields.res = {'val', 'width', 'nSquares'};
circuit.allowedPrintFields.cap = {'val'};
circuit.allowedPrintFields.ind = {'val'};
circuit.allowedPrintFields.v = {'val'};
circuit.allowedPrintFields.i = {'val'};

if exist('nmosTable', 'var') & isstruct(nmosTable)
  circuit = cirInitMosDefaults(circuit, nmosTable, 'nmos', nmosTableName);
  % this is an inline function (see below)
end
if exist('pmosTable', 'var') & isstruct(pmosTable)
  circuit = cirInitMosDefaults(circuit, pmosTable, 'pmos', pmosTableName);
  % this is an inline function (see below)
end


circuit.defaults.res.fieldNames = {'val', 'paramFields', 'paramPrintNames'};
circuit.defaults.res.fieldValues{1} = NaN;

circuit.defaults.cap.fieldNames = {'val', 'paramFields', 'paramPrintNames'};
circuit.defaults.cap.fieldValues{1} = NaN;

circuit.defaults.ind.fieldNames = {'val', 'paramFields', 'paramPrintNames'};
circuit.defaults.ind.fieldValues{1} = NaN;

circuit.defaults.v.fieldNames = {'val', 'paramFields', 'paramPrintNames'};
circuit.defaults.v.fieldValues{1} = NaN;

circuit.defaults.i.fieldNames = {'val', 'paramFields', 'paramPrintNames'};
circuit.defaults.i.fieldValues{1} = NaN;

if exist('designkitName', 'var')
  if exist('simulator', 'var')
    designkit = cirDesignkitInit(designkitName, simulator);
  else
    designkit = cirDesignkitInit(designkitName, 'none');
  end
  circuit = cirInitParamPrintNames(circuit, designkit.printParamTuple);
  circuit.designkit = designkit;
  % check on consistency between designkit and tables:
  if exist('nmosTable', 'var') & isstruct(nmosTable)
    cirDesignkitAndTableCheck(circuit.designkit, nmosTable); 
    % this is an inline function, see below 
  end
  if exist('pmosTable', 'var') & isstruct(pmosTable)
    cirDesignkitAndTableCheck(circuit.designkit, pmosTable); 
    % this is an inline function, see below 
  end
  
end


% initially, the circuit does not contain any element:
circuit.nElements = 0;
% initialization of number of nmos transistors, number of pmos transistors,
% etc. to zero:
eltypes = cirSupportedEltypes;
for i = 1:length(eltypes)
  circuit.(strcat('n', eltypes{i})) = 0;
end

circuit = cirElementsCheckIn(circuit, elementList, hierarchy);

% END OF FUNCTION CIRINIT

% definition of inline functions:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function circuit = cirInitMosDefaults(circuit, table, mostype, tableName)
% mostype is either the string 'nmos' or the string 'pmos'
circuit.defaults.(mostype).fieldNames = {'table', 'lg', 'geo', ...
    'paramFields', 'paramPrintNames', 'lsos', 'lsod', 'lsogs', 'lsogd', ...
    'mult'};
% we now specify the default values for the fields defined in the list
% circuit.defaults.(mostype).fieldNames: 

% default value for table:
circuit.defaults.(mostype).fieldValues{1} = tableName;

% default length:
circuit.defaults.(mostype).fieldValues{2} = tableLmin(table);

% default value for geo: 1
% this means even number of fingers, more source fingers than drain fingers
circuit.defaults.(mostype).fieldValues{3} = 1;

% default value for lsos:
if (isModelParam('wco', table) & isModelParam('eactco', table) & ...
      isModelParam('dcopss', table))
  circuit.defaults.(mostype).fieldValues{6} = mosLsoS(tableModelParam('wco', ...
      table), tableModelParam('eactco', table), tableModelParam('dcopss', table));
elseif (isModelParam('hdif', table))
  circuit.defaults.(mostype).fieldValues{6} = 2*tableModelParam('hdif', ...
      table);
else
  error('Too few parameters available for computation of lsos');
end

% default value for lsod:
if (isModelParam('wco', table) & isModelParam('eactco', table) & ...
      isModelParam('dcopsd', table))
  circuit.defaults.(mostype).fieldValues{7} = mosLsoD(tableModelParam('wco', ...
      table), tableModelParam('eactco', table), tableModelParam('dcopsd', table));
elseif (isModelParam('hdif', table))
  circuit.defaults.(mostype).fieldValues{7} = 2*tableModelParam('hdif', ...
      table);
else
  error('Too few parameters available for computation of lsod');
end

% default value for lsogs:
if (isModelParam('wco', table) & isModelParam('dcopss', table))
  circuit.defaults.(mostype).fieldValues{8} = mosLsoGs(tableModelParam('wco', ...
      table),  tableModelParam('dcopss', table));
elseif (isModelParam('hdif', table))
  circuit.defaults.(mostype).fieldValues{8} = 2*tableModelParam('hdif', ...
      table);
else
  error('Too few parameters available for computation of lsogs');
end

% default value for lsogd:
if (isModelParam('wco', table) & isModelParam('dcopsd', table))
  circuit.defaults.(mostype).fieldValues{9} = ...
      mosLsoGd(tableModelParam('wco', table), tableModelParam('dcopsd', table));
elseif (isModelParam('hdif', table))
  circuit.defaults.(mostype).fieldValues{9} = 2*tableModelParam('hdif', ...
      table);
else
  error('Too few parameters available for computation of lsogd');
end

% default value for mult:
circuit.defaults.(mostype).fieldValues{10} = 1;

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function circuit = cirInitParamPrintNames(circuit, tuple)
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

if length(tuple) ~= 3*length(cirSupportedEltypes) - 3
  error('wrong length of the tuple (value %d)', length(tuple));
end

for i = 1:3:length(tuple)-2
  if not(ischar(tuple{i}))
    error('circuit element type is not a string');
  end
  if not(ismember(tuple{i}, cirSupportedEltypes))
    error('unknown circuit element type %s', tuple{i})
  end
  if strcmp(tuple{i}, 'subckt')
    errorstring = strcat('Initialization of print parameters for subcircuits', ...
	'is not done in this way');
    error(errorstring);
  end
end

for i = 2:3:length(tuple)-1
  if length(tuple{i}) ~= length(tuple{i+1})
    error('length of element %d is not equal to the length of element %d\n', ...
	i, i+1);
  end
end

for i = 1:3:length(tuple)
  eltype = tuple{i};
  for j=1:length(tuple{i+1})
  if not(ismember(tuple{i+1}{j}, cirAllowedPrintFields(circuit, eltype)))
    error('item %s is not an allowed print field for elements of type %s', ...
	tuple{i+1}{j}, eltype); 
  end
  switch eltype
  case {'nmos', 'pmos'}
    circuit.defaults.(eltype).fieldValues{4} = tuple{i+1};
    circuit.defaults.(eltype).fieldValues{5} = tuple{i+2};
  case {'res', 'cap', 'ind', 'v', 'i'}
    circuit.defaults.(eltype).fieldValues{2} = tuple{i+1};
    circuit.defaults.(eltype).fieldValues{3} = tuple{i+2};
  otherwise    
    error('Unknown element type');
end

end
end

function cirDesignkitAndTableCheck(designkit, table)
if designkit.tech.node ~= tableTechNode(table)
  fprintf(1, 'Warning: The technology node of the design kit is %sm\n', ...
      eng(designkit.tech.node));
  fprintf(1, 'This is different from the technology node of the table'); 
  fprintf(1, '\nof the %s-MOS transistor which is %sm\n\n', ...
      tableType(table), eng(tableTechNode(table)));
end
if not(strcmpi(designkit.tech.source, tableTechSource(table)))
  fprintf(1, 'Warning: the technology source of the design kit is %s\n', ...
      designkit.tech.source);
  fprintf(1, 'This is different from the source of the table');
  fprintf(1, '\nof the %s-MOS transistor which is %s\n\n', ...
      tableType(table), tableTechSource(table));
end


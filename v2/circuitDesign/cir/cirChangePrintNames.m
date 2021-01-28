function circuit = cirChangePrintNames(circuit, eltype, nameList, printNameList)
%CIRCHANGEPRINTNAMES adapts the list of fields that will be written into a 
%    circuit simulation file.
%
%    cirChangePrintNames(CIRCUIT, ELTYPE, NAMELIST, PRINTNAMELIST) modifies two
%    cell arrays:
%    1. the one that contains the parameters of a circuit element of a given
%    type ELTYPE (must be a string) is set equal to the argument NAMELIST;
%    2. the one that contains for each of these parameters the string that is
%    written into the simulation file, is set equal to the argument
%    PRINTNAMELIST. paramFields and the paramPrintNames, which
%    are strings that are the fieldnames of a circuit element of type ELTYPE
%    when this circuit element is defined with cirElementDef.
%
%    This function changes the fields circuit.defaults.(eltype).fieldValues{5}
%    and circuit.defaults.(eltype).fieldValues{6}. If this function is called,
%    it must be done before any circuit element is defined.
%
%    See also cirInit, cirElementDef, cirDefaultValues. 
%
%    EXAMPLE:
%    
%   cir = cirChangePrintNames(cir, 'nmos', {'w', 'lg', 'nFingers'}, ...
%     {'Wnom', 'Lnom', 'FOLD'})
%
%   With this example, we determine for circuit "cir" that in the circuit simulation file which
%   will be written or modified, the width (w), the length (lg) and the number
%   of fingers of an n-MOS will be written/modified. For example, in a circuit
%   netlist (Spectre format) that contains n-MOS transistor Mn1, for which we
%   have computed in Matlab w=10e-6, lg=0.2e-6 and nFingers = 5, the above
%   example would yield a line in the circuit simulation file of the form
%
%   Mn1 (drain gate source bulk) modelname Wnom=10e-6 Lnom=0.2e-6 FOLD=5
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%


cirCheckEltype(eltype);
if length(nameList) ~= length(printNameList)
  error('lengths of NAMELIST and PRINTNAMELIST are unequal');
end

for i = 1:length(nameList)
  if not(ismember(nameList{i}, cirAllowedPrintFields(circuit, eltype)))
    error('item %s is not an allowed print field', nameList{i});
  end
end

switch eltype
  case {'nmos', 'pmos'}
    circuit.defaults.(eltype).fieldValues{5} = nameList;
    circuit.defaults.(eltype).fieldValues{6} = printNameList;
  case {'res', 'cap', 'ind', 'v', 'i'}
    circuit.defaults.(eltype).fieldValues{2} = nameList;
    circuit.defaults.(eltype).fieldValues{3} = printNameList;
  case 'subckt'
    circuit.
  otherwise    
    error('Unknown element type');
end


    

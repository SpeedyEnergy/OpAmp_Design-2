function circuit = cirElementsCheckIn(circuit, elementList, top)
%Definition of circuit elements for a given circuit.
%
%  CIRELEMENTSCHECKIN(CIRCUIT, ELEMENTLIST, TOP) defines in the given CIRCUIT
%  for a circuit element with name '<name>', a field CIRCUIT.<name>. This field
%  is initialized as a circuit element of a certain type, using the function
%  cirElementDef. The following types of circuit elements are allowed:
%  voltage sources (type "v"), current sources (type "i"), n-MOS transistors
%  (type "nmos"), p-MOS transistors (type "pmos"), capacitors (type "cap"),
%  resistors (type "res"), inductors (type "ind") and subcircuits (type
%  "subckt"). The circuit elements must be specified with the argument
%  ELEMENTLIST, which is a structure with the following fields:
%  - if independent voltage sources must be defined, then ELEMENTLIST should
%    have a field v. If so, then ELEMENTLIST.v should be a cell array, that
%    contains the names of all independent voltage sources. If no independent
%    voltage sources have to be defined, then either the field ELEMENTLIST.v is
%    not defined, or it is equal to the empty cell array {}; 
%  - if independent current sources must be defined, then ELEMENTLIST should
%    have a field i. If so, then ELEMENTLIST.i should be a cell array, that
%    contains the names of all independent current sources. If no independent
%    current sources have to be defined, then either the field ELEMENTLIST.i is
%    not defined, or it is equal to the empty cell array {};
%  - if n-MOS transistors must be defined, then ELEMENTLIST should have a field
%    nmos. If so, then ELEMENTLIST.nmos should be a cell array, that
%    contains the names of all n-MOS transistors. If no n-MOS transistors have
%    to be defined, then either the field ELEMENTLIST.nmos is not defined, or
%    it is equal to the empty cell array {}; 
%  - idem for the circuit elements of type "pmos", "cap", "res", "ind" and
%    "subckt".
%
%  If the function CIRELEMENTSCHECKIN is called for the top circuit (which is
%  always the case for a circuit without hierarchy), then the argument TOP must
%  have a nonzero value. In this case, in addition to defining circuit elements
%  of the form CIRCUIT.<name>, the variable <name> is also copied to the base
%  workspace. This has been provided for easier manipulation of the circuit
%  elements without too much typing overhead and for a better readibility of
%  the code. For example for a MOS transistor with name 'Mn1' in the circuit
%  with name 'opamp', manipulating the variable Mn1 in the base workspace, is
%  more elegant than having to manipulate opamp.Mn1.
%  Important remark: if the variable <name> is used for the sizing of the
%  circuit, then, after the sizing of <name>, all fields of <name> should be
%  copied to CIRCUIT.<name>. This is done with the function
%  cirElementsCheckOut. 
%
%  See also cirElementDef, cirElementsCheckOut
%
%  EXAMPLES: 
%
%  Let's define an elementary OTA. The name of the corresponding variable will
%  be "elOta". First we define a list of circuit elements, grouped per element
%  type: 
%     elementList.v = {'Vdd', 'VinN', 'VinP', 'Vbias'};
%     elementList.nmos = {'Mn1a', 'Mn1b', 'Mnbias'};
%     elementList.pmos = {'Mp1a', 'Mp1b'};
%     elementList.i = {};

%  Then the following command
%
%     elOta = cirElementsCheckIn(elOta, elementList, 'top')
%
%  defines the following fields
%
%     elOta.Vdd, elOta.VinN, elOta.VinP, elOta.Vbias, elOta.Mn1a, elOta.Mn1b,
%     elOta.Mnbias, elOta.Mp1a and elOta.Mp1b.
%
%  These fields are initialized, according to the element type, using the
%  function "cirEelementDef".
%  In addition, the following variables are defined in the base workspace:
%     Vdd, VinN, VinP, Mn1a, Mn1b, Mp1a and Mp1b.
%  These variables are copies of the variables elOta.Vdd, elOta.VinN, ...
%  The variables Vdd, VinN, VinP, Mn1a, Mn1b, Mp1a and Mp1b can now be used for
%  sizing of the circuit. After sizing, these variables get new fields and
%  existing fields can have new values. After the sizing of the circuit, the
%  variables Vdd, VinN, VinP, Mn1a, Mn1b, Mp1a and Mp1b can be copied to the
%  variables elOta.Vdd, elOta.VinN, elOta.VinP, elOta.Mn1a, ... using the
%  function cirElementsCheckOut. 
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

debug = 0;
if debug
  disp(elementList)
end

if isfield(elementList, 'v')
  for i = 1:length(elementList.v)
    circuit = cirElementDef(circuit, 'v', elementList.v{i});
  end
end

if isfield(elementList, 'i')
  for i = 1:length(elementList.i)
    circuit = cirElementDef(circuit, 'i', elementList.i{i});
  end
end


if isfield(elementList, 'nmos')
  for i = 1:length(elementList.nmos)
    circuit = cirElementDef(circuit, 'nmos', elementList.nmos{i});
  end
end

if isfield(elementList, 'pmos')
  for i = 1:length(elementList.pmos)
    circuit = cirElementDef(circuit, 'pmos', elementList.pmos{i});
  end
end

if isfield(elementList, 'cap')
  for i = 1:length(elementList.cap)
    circuit = cirElementDef(circuit, 'cap', elementList.cap{i});
  end
end

if isfield(elementList, 'res')
  for i = 1:length(elementList.res)
    circuit = cirElementDef(circuit, 'res', elementList.res{i});
  end
end

if isfield(elementList, 'ind')
  for i = 1:length(elementList.ind)
    circuit = cirElementDef(circuit, 'ind', elementList.ind{i});
  end
end

if isfield(elementList, 'subckt')
  for i = 1:length(elementList.subckt)
    circuit = cirElementDef(circuit, 'subckt', elementList.subckt{i});
  end
end

if debug
  fprintf(1, 'circuit is\n');
  disp(circuit);
end

if top
  for i = 1:circuit.nv
    assignin('base', elementList.v{i}, circuit.(circuit.vList{i}));
  end
  for i = 1:circuit.ni
    assignin('base', elementList.i{i}, circuit.(circuit.iList{i}));
  end
  for i = 1:circuit.nnmos
    assignin('base', elementList.nmos{i}, circuit.(circuit.nmosList{i}));
  end
  for i = 1:circuit.npmos
    assignin('base', elementList.pmos{i}, circuit.(circuit.pmosList{i}));
  end
  for i = 1:circuit.ncap
    assignin('base', elementList.cap{i}, circuit.(circuit.capList{i}));
  end
  for i = 1:circuit.nres
    assignin('base', elementList.res{i}, circuit.(circuit.resList{i}));
  end
  for i = 1:circuit.nind
    assignin('base', elementList.ind{i}, circuit.(circuit.indList{i}));
  end
  for i = 1:circuit.nsubckt
    assignin('base', elementList.subckt{i}, circuit.(circuit.subcktList{i}));
  end

end

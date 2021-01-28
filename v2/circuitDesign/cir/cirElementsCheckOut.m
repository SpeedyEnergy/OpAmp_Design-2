function circuit = cirElementsCheckOut(circuit)
%CIRELEMENTSCHECKOUT copies circuit element with name <x> to <circuit>.<x>
%
% See also cirElementsCheckIn
%
%  EXAMPLES:
%
% we have a circuit with name 'commSource'. The circuit elements are
% commSource.Mn1, commSource.Rl, commSource.Vdd and commSource.Vin
% With the function cirElementsCheckIn these elements have been copied to the
% variables Mn1, Rl, Vdd and Vin. In the Matlab code that sizes these elements,
% the fields of these elements get meaningful values, e.g. Vin.val gets a
% value, ...
% After this sizing we can use cirElementsCheckOut to copy the fields of Mn1,
% Rl, Vdd and Vin back to commSource.Mn1, commSource.Rl, commSource.Vdd and
% commSource.Vin. This is done by the command
%
% commSource = cirElementsCheckOut(commSource)
%
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

for i = 1:circuit.nnmos
  circuit.(circuit.nmosList{i}) =  evalin('base', circuit.nmosList{i});
end
for i = 1:circuit.npmos
  circuit.(circuit.pmosList{i}) =  evalin('base', circuit.pmosList{i});
end
for i = 1:circuit.nv
  circuit.(circuit.vList{i}) =  evalin('base', circuit.vList{i});
end
for i = 1:circuit.ni
  circuit.(circuit.iList{i}) =  evalin('base', circuit.iList{i});
end
for i = 1:circuit.ncap
  circuit.(circuit.capList{i}) =  evalin('base', circuit.capList{i});
end
for i = 1:circuit.nres
  circuit.(circuit.resList{i}) =  evalin('base', circuit.resList{i});
end
for i = 1:circuit.nind
  circuit.(circuit.indList{i}) =  evalin('base', circuit.indList{i});
end
for i = 1:circuit.nsubckt
  circuit.(circuit.subcktList{i}) =  evalin('base', circuit.subcktList{i});
end

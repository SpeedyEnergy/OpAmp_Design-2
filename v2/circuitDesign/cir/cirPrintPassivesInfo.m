function cirPrintPassivesInfo(fid, circuit)
% CIRPRINTPASSIVESINFO print element values and other relevant
% information of passive components in a circuit

fprintf(fid, '\n******************************************\n');
fprintf(fid, 'circuit %s   passive components \n', circuit.name); 
fprintf(fid, '*******************************************\n');

if circuit.nres > 0 
  % print resistor info:
  fprintf(fid, 'Resistors:\n');
  for i = 1:circuit.nres
    cirPrintResistor(fid, circuit.(circuit.resList{i}));
  end
  fprintf(fid, '\n');
end

if circuit.ncap > 0 
  % print capacitor info:
  fprintf(fid, 'Capacitors:\n');
  for i = 1:circuit.ncap
    cirPrintCapacitor(fid, circuit.(circuit.capList{i}));
  end
  fprintf(fid, '\n');
end

if circuit.nind > 0 
  % print inductor info:
  fprintf(fid, 'Inductors:\n');
  for i = 1:circuit.nind
    cirPrintInductor(fid, circuit.(circuit.indList{i}));
  end
  fprintf(fid, '\n');
end

fprintf(fid, '*****************************************\n\n');

function cirPrintResistor(fid, element)

fprintf(fid, '%s: value = %sOhm', element.name, eng(element.val));
if isfield(element, 'width')
  fprintf(fid, '\twidth = %sm', eng(element.width));
end
if isfield(element, 'nSquares')
  fprintf(fid, '\tnSquares = %d', element.nSquares);
end

fprintf(fid, '\n');

% end of function cirPrintResistor

function cirPrintCapacitor(fid, element)

fprintf(fid, '%s: value = %sF', element.name, eng(element.val));

fprintf(fid, '\n');

% end of function cirPrintCapacitor

function cirPrintInductor(fid, element)

fprintf(fid, '%s: value = %sH', element.name, eng(element.val));

fprintf(fid, '\n');

% end of function cirPrintInductor


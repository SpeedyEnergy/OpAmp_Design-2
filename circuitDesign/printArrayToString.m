function result = printArrayToString(array, initPosition, lineLength, ...
    interString, blanks)
%PRINTARRAYTOSTRING makes a string including carriage that includes a printout
%    of an array
%
%   The array must consist of numbers only or strings only.
if isnumeric(array(1))
  format = '%g';
else
  format = '%s';
end

interStringLength = length(interString);
position = initPosition;
result = '';
if (position == 0)
  for j=1:blanks
    result = sprintf(' %s', result);
  end
  position = blanks;
end

if iscell(array)
  printString = sprintf(format, char(array(1)));
else
  printString = sprintf(format, array(1));
end
lengthPrintString = length(printString);
result = sprintf('%s%s%s', result, printString, interString);
position = position + lengthPrintString + interStringLength;

for i = 2:length(array)-1
  if position > lineLength
    result = sprintf('%s\n', result);
    for j=1:blanks
      result = sprintf('%s ', result);
    end
    position = blanks;
  end

  if iscell(array)
    printString = sprintf(format, char(array(i)));
  else
    printString = sprintf(format, array(i));
  end
  lengthPrintString = length(printString);
  result = sprintf('%s%s%s', result, printString, interString);
  position = position + lengthPrintString + interStringLength;
end

if position > lineLength
  result = sprintf('%s\n', result)
end

if iscell(array)
  printString = sprintf(format, char(array(length(array))));
else
  printString = sprintf(format, array(length(array)));
end

result = sprintf('%s%s', result, printString);


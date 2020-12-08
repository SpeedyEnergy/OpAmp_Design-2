function resultIndex = binSearch(value, array)
%BINSEARCH binary search of the index in a monotonic array that corresponds 
% to the element which is closest to the given value 
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%
finalIndex = length(array);
beginIndex = 1;
i = 0;
while i < 100
  medianIndex = beginIndex + fix((finalIndex - beginIndex)/2);
  if (finalIndex - beginIndex == 1) | (value == array(medianIndex))
    resultIndex = medianIndex;
    return;
  elseif value < array(medianIndex)
    finalIndex = medianIndex;
  else
    beginIndex = medianIndex;
  end
  i = i+1;
end
if i == 100
  error('we have gone through 100 searches in this binary search...');
end

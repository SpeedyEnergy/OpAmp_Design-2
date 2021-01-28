function result = noiseSourceIsFlicker(contrib)
if strcmp(contrib.source.type, 'flicker')
  result = 1;
else
  result = 0;
end

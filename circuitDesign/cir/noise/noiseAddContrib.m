function contrib = noiseContrib(name, value, type, outName, tfDc, zeros, ...
    poles, eqInName)
%NOISECONTRIB constructor function for a noise contribution

contrib.source.name = name;
contrib.source.val = value;
switch type
  case 'flicker'
    contrib.source.type = 1;
  case 'thermal'
    contrib.source.type = 0;
  otherwise
    error('unknown type of noise source');
end

contrib.tf.dcVal = tfDc;
if iscell(poles)
  contrib.tf.poles = poles;
else
  error('poles should be specified as a cell array');
end

if iscell(zeros)
  contrib.tf.zeros = zeros;
else
  error('zeros should be specified as a cell array');
end

contrib.out.name = outName;
contrib.eqIn.name = eqInName;



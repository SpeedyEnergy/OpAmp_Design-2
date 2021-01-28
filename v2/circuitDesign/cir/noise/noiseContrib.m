function contrib = noiseContrib(name, value, type, outName, outType, tfDc, ...
    zeros, poles)
%NOISECONTRIB constructor function for a noise contribution

contrib.source.name = name;
contrib.source.val = value;
switch type
  case 'flicker'
    contrib.source.type = 'flicker';
  case 'thermal'
    contrib.source.type = 'thermal';
  otherwise
    error('unknown type of noise source');
end

if not(strcmp(outType, 'v')) & not(strcmp(outType, 'i'))
  error('Wrong type of the output quantity, namely %s', outType)
end

contrib.tf.dcVal = tfDc;
contrib.tf.poles = poles;
contrib.tf.zeros = zeros;

contrib.out.name = outName;
contrib.out.type = outType;



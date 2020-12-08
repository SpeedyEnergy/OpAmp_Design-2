function contrib = noiseEnterSource(name, value, type, outName, outType, tf)
%NOISEENTERSOURCE constructor function for a noise source and its contribution 
%    to the total noise 
%
%   CONTRIB = NOISEENTERSOURCE(NAME, VALUE, TYPE, OUTNAME, OUTTYPE, TF) returns
%   a structure CONTRIB, that corresponds to data of one single noise
%   source. The structure CONTRIB has the following fields:
%   1. source
%      This has the following fields:
%       - name
%         This is a string, representing the name of the source (e.g. 'Mp1
%         thermal')
%       - val
%         This is the value of the noise source. This can be a tabulated value
%         (e.g. Mp1.di2_id) or an equivalent input noise voltage (then one
%         should use Mp1.di2_id / Mp1.gm^2)
%       - type
%         This is a string, which is either 'flicker' or 'thermal'
%   2. tf
%      This is the transfer function from the noise source to the circuit
%      output. This circuit output is specified by its name, which is the
%      argument OUTNAME (a string, e.g. 'vout'), and its type OUTTYPE. The type
%      is either 'v' (when the output of interest is a voltage) or 'i' (when
%      the output of interest is a current).
%   3. out
%      This has two fields initially:
%      - name
%        This is the name of the circuit output (= OUTNAME, see above).
%      - type
%        This is the type of the circuit output ('v' or 'i', see above)
%   When the noise power spectral density or noise figure is computed (see
%   function noiseCompute), then extra fields are added to CONTRIB.
%   - the field out gets as extra field "psd"
%   - if the equivalent input noise power spectral density is computed with
%     noiseCompute, then the field eqIn is added. See noiseCompute for more
%     details.
%   - if the spot noise figure is computed with noiseCompute, then the field
%   "nfLin" is added. See noiseCompute for more details.

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

contrib.tf = tf;

contrib.out.name = outName;
contrib.out.type = outType;



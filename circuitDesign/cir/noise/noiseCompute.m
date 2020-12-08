function noise = noiseCompute(contribArray, outName, outType, tf, freq, varargin)
%NOISECOMPUTE computation of the output noise spectral density and other 
%   derived quantities.
%
%  NOISE = noiseCompute(CONTRIBARRAY, OUTNAME, OUTTYPE, TF, FREQ) computes the
%  power spectral density at the output of a circuit. This output is specified
%  with a name (string OUTNAME), a type (string OUTTYPE) and a transfer
%  function TF. This transfer function is a structure with the following three
%  fields:  
%    - dcVal: the value of the transfer function for s = jw = 0
%    - poles: an array containing the value of the poles that should be
%             taken into account.
%    - zeros: an array containing the value of the zeros that should be
%             taken into account.
%  The total output power spectral density depends on different
%  contributions, which are given in the array CONTRIBARRAY. Each contribution
%  is a structure, with the datastructure that has been constructed using the
%  function noiseContrib.
%  The units of the resulting power spectral density are
%  A^2/Hz or V^2/Hz, depending on the type of the output of interest. The power
%  spectral density is computed as a function of frequency, which is specified
%  by the argument FREQ. This is an array of frequency points.
%  The results are returned as part of the structure NOISE.
%
%  NOISE = noiseCompute(CONTRIBARRAY, OUTNAME, OUTTYPE, TF, FREQ, EQINNAME,
%  EQINTYPE) computes in addition the equivalent input noise power spectral
%  density at the input, which is specified by a string (input argument
%  EQINNAME) and a type (string EQINTYPE), which is either 'v' or 'i'.
%
%  NOISE = noiseCompute(CONTRIBARRAY, OUTNAME, OUTTYPE, TF, FREQ, EQINNAME,
%  EQINTYPE, SOURCEINDEX) computes in addition the spot noise figure as a
%  function of frequency. Here, SOURCEINDEX is the index in the array
%  CONTRIBARRAY that corresponds to the noise source of the source impedance.
%
%  The datastructure of the return value NOISE is a structure with the
%  following fields:
%
%  1. contrib
%     This is an array of contributions: contrib(1), contrib(2),
%     ... contrib(i), ... It is a copy of the argument CONTRIBARRAY, with, for
%     each contribution contrib(i), the addition of the following fields: 
%     - contrib(i).out.psd
%       this is the contribution of the i-th noise source to the power spectral
%       density of the noise at the output;
%     - contrib(i).eqIn.psd
%       this is the contribution of the i-th noise source to the power spectral
%       density of the noise at the output, referred back to the input;
%     - contrib(i).nfLin
%       this is the contribution of the i-th noise source to the noise figure
%       (linear scale, not in dB);
%  2. freq
%     This is an array of frequency points, describing the frequency band over
%     which noise power spectral density and integrated noise are computed
%  3. nContribs
%     number of contributions
%  4. out
%     This is in turn a structure with the following fields:
%      - name
%        this is the name (= string) of the output quantity;
%      - type
%        this is the type (= string) of the output quantity. This is either 'v'
%        (for a voltage) or 'i' (for a current);
%      - psd
%        this is an array with the same length as the field "freq", with, for
%        every corresponding frequency point, the output noise power spectral
%        density;
%  5. eqIn
%     This is a structure that is constructed if equivalent input noise is to
%     be computed. The structure has the following fields:
%      - name
%        this is the name (= string) of the input signal
%      - type
%        this is the type (= string) of the input signal. This is either 'v'
%        (for a voltage) or 'i' (for a current);
%      - psd
%        this is an array with the same length as the field "freq", with, for
%        every corresponding frequency point, the input-referred noise power
%        spectral density
%      - rms
%        square root of the input-referred noise power, which is the integral
%        of the input-referred noise power spectral density over the frequency
%        band of interest, which is specified with "freq"
%  6. tf
%     This is the transfer function, which is a structure with three fields:
%      - dcVal: the value of the transfer function for s = jw = 0
%      - poles: an array containing the value of the poles that should be
%             taken into account.
%      - zeros: an array containing the value of the zeros that should be
%             taken into account.
%  7. nfLin
%     This field is computed when the noise figure is to be computed. It is an
%     array with the same length as the field "freq", with for every
%     corresponding frequency point the noise figure (linear scale)
%
%  Added by Jonathan:
%     The following fields are added to the noise structure:
%         - For the contributors 'averagedOut' (V^2/Hz), 'integratedOut' (V^2),
%         'averagedEqIn' (V^2) and 'integratedEqIn' (V^2/Hz) are resp. the
%         averaged output noise contribution over the bandwidth determined in
%         freq, the total integrated noise, the averaged input referred
%         noise and the integrated input referred noise.
%         - For the general structure similar parameters are added, being
%         the total noise values.
%
%  (c) IMEC, 2005
%  IMEC confidential 
%


debug = 0;
noise.contrib = contribArray;
noise.freq = freq;
noise.tf = tf;
noise.out.name = outName;
noise.nContribs = length(contribArray);
if not(strcmp(outType, 'v')) & not(strcmp(outType, 'i'))
  error('Wrong type of the output quantity, namely %s', outType)
end
noise.out.type = outType;

for i = 1:noise.nContribs
  if not(strcmp(noise.contrib(i).out.name, outName))
    errorString = strcat(...
	sprintf('output of noise source %d is %s.\n', i, ...
	noise.contrib(i).out.name), ...
	sprintf('This differs from circuit output %s\n', outName));
    error(errorString);
  end
  if not(strcmp(noise.contrib(i).out.type, noise.out.type))  
    errorString = strcat(...
	sprintf('type output of noise source %d is %s.\n', i, ...
	noise.contrib(i).out.type), ...
	sprintf('This differs from type of circuit output %s\n', noise.out.type));
    error(errorString);
  end
end


if nargin > 5
  if nargin == 6
    error('The input signal is specified but not its type');
  end
  noise.eqIn.name = varargin{1};
  noise.eqIn.type = varargin{2};
  computeEqIn = 1;
else
  computeEqIn = 0;
end

if nargin > 7
  computeNf = 1;
  sourceIndex = varargin{3};
  if not(isnumeric(sourceIndex)) | (sourceIndex < 1) | (sourceIndex > ...
	noise.nContribs)
    error('Wrong index for the noise of the source impedance');
  end
else
  computeNf = 0;
end

for i = 1:noise.nContribs
  % in case of tf values:
  if isfield(noise.contrib(i).tf, 'values')
    noise.contrib(i).out.psd = (abs(noise.contrib(i).tf.values)).^2;
  % in case of poles and zeros:
  elseif isfield(noise.contrib(i).tf, 'poles') & isfield(noise.contrib(i).tf, ...
	'zeros') 
    noise.contrib(i).out.psd =  noise.contrib(i).source.val * ...
	noise.contrib(i).tf.dcVal^2 * ones(size(freq));
    if debug
      fprintf(1, '\ncontribution of %s: initially it has a length of %d\n', ...
	  noise.contrib(i).source.name)
    end
    % handle poles
    for j = 1:length(noise.contrib(i).tf.poles)
      a = real(noise.contrib(i).tf.poles(j).val);
      b = imag(noise.contrib(i).tf.poles(j).val);
      if debug
	fprintf(1, ...
	    'pole %d for output noise: real part = %g, imag. part = %g\n', ...
	    j, a, b);
      end
      noise.contrib(i).out.psd = (a^2 + b^2) * noise.contrib(i).out.psd ./ (a^2 ...
	  + (b + 2*pi*freq).^2);
    end
    % handle zeros
    for j = 1:length(noise.contrib(i).tf.zeros)
      a = real(noise.contrib(i).tf.zeros(j).val);
      b = imag(noise.contrib(i).tf.zeros(j).val);
      if debug
	fprintf(1, ...
	    'zero %d for output noise: real part = %g, imag. part = %g\n', ...
	    j, a, b);
      end
      noise.contrib(i).out.psd = 1/(a^2 + b^2) * noise.contrib(i).out.psd .* ...
	  (a^2 + (b + 2*pi*freq).^2);
    end
  % in case of numerator and denumerator:
  elseif isfield(noise.contrib(i).tf, 'denomCoeffs') & ...
	isfield(noise.contrib(i).tf, 'numerCoeffs') 
    noise.contrib(i).out.psd = (abs(polyval(noise.contrib(i).tf.numerCoeffs, ...
	2*pi*sqrt(-1)*freq))).^2 ./ ...
	(abs(polyval(noise.contrib(i).tf.denomCoeffs, 2*pi*sqrt(-1)*freq))).^2;
	
  end
  
  if noiseSourceIsFlicker(noise.contrib(i))
    fprintf(1, 'noise source %s is a 1/f noise source\n', ...
	noise.contrib(i).source.name);
    noise.contrib(i).out.psd = noise.contrib(i).out.psd ./ freq;
  end

  if computeEqIn
    if isfield(noise.tf, 'values')
      noise.contrib(i).eqIn.psd = noise.contrib(i).out.psd ./ ...
	  (abs(noise.tf.values)).^2;
    elseif isfield(noise.tf, 'zeros') & isfield(noise.tf, 'poles')
      noise.contrib(i).eqIn.psd = noise.contrib(i).out.psd ./ noise.tf.dcVal^2;
      for j = 1:length(noise.tf.zeros)
	a = real(noise.tf.zeros(j).val);
	b = imag(noise.tf.zeros(j).val);
	if debug
	  fprintf(1, ...
	      'zero %d for eq. input noise: real part = %g, imag. part = %g\n', ...
	      j, a, b);
	end
	noise.contrib(i).eqIn.psd = (a^2 + b^2) * noise.contrib(i).eqIn.psd ./ ...
	    (a^2 + (b + 2*pi*freq).^2);
      end
      for j = 1:length(noise.tf.poles)
	a = real(noise.tf.poles(j).val);
	b = imag(noise.tf.poles(j).val);
	if debug
	  fprintf(1, ...
	      'pole %d for eq. input noise: real part = %g, imag. part = %g\n', ...
	      j, a, b);
	end
	noise.contrib(i).eqIn.psd = 1/(a^2 + b^2) * noise.contrib(i).eqIn.psd .* ...
	    (a^2 + (b + 2*pi*freq).^2);
      end
    elseif isfield(noise.tf, 'numerCoeffs') & isfield(noise.tf, ...
	  'denomCoeffs')
      noise.contrib(i).eqIn.psd = noise.contrib(i).out.psd ./ ...
	  (abs(polyval(noise.tf.numerCoeffs, 2*pi*sqrt(-1)*freq))).^2 .* ...
	  (abs(polyval(noise.tf.denomCoeffs, 2*pi*sqrt(-1)*freq))).^2
    end
  end
end


% computation of the total output noise power spectral density:
noise.out.psd = zeros(size(freq));
if debug
  fprintf(1, 'fields of noise.out:\n');
  disp(fieldnames(noise.out));
end
for i = 1:noise.nContribs
  noise.out.psd = noise.out.psd + noise.contrib(i).out.psd;
  if debug
    fprintf(1, 'fields of noise.contrib(%d).out:\n', i);
    disp(fieldnames(noise.contrib(i).out));
  end
end

% computation of the total input-referred noise power spectral density:
if computeEqIn
  noise.eqIn.psd = zeros(size(freq));
  for i = 1:noise.nContribs
    noise.eqIn.psd = noise.eqIn.psd + noise.contrib(i).eqIn.psd;
  end
end

% computation of the noise figure:
if computeNf
  noise.nfLin = 0;
  for i = 1:noise.nContribs
    noise.contrib(i).nfLin = noise.contrib(i).out.psd ./ ...
	noise.contrib(sourceIndex).out.psd;
    noise.nfLin = noise.nfLin + noise.contrib(i).nfLin;
  end
end

% computation of the integrated and averaged noise:
% this section only works for poles and zeros!
nPoints = 10000 ;
maxFreq = max(freq) ;
minFreq = min(freq) ;
freqBand = maxFreq - minFreq ;
stepSize = freqBand / nPoints ;
freqLin = minFreq:stepSize:maxFreq ;
for i = 1:noise.nContribs
  % in case of poles and zeros:
  if isfield(noise.contrib(i).tf, 'poles') & isfield(noise.contrib(i).tf, ...
	'zeros') 
    noise.contrib(i).linout.psd =  noise.contrib(i).source.val * ...
	noise.contrib(i).tf.dcVal^2 * ones(size(freqLin));
    % handle poles
    for j = 1:length(noise.contrib(i).tf.poles)
      a = real(noise.contrib(i).tf.poles(j).val);
      b = imag(noise.contrib(i).tf.poles(j).val);
      noise.contrib(i).linout.psd = (a^2 + b^2) * noise.contrib(i).linout.psd ./ (a^2 ...
	  + (b + 2*pi*freqLin).^2);
    end
    % handle zeros
    for j = 1:length(noise.contrib(i).tf.zeros)
      a = real(noise.contrib(i).tf.zeros(j).val);
      b = imag(noise.contrib(i).tf.zeros(j).val);
      noise.contrib(i).linout.psd = 1/(a^2 + b^2) * noise.contrib(i).linout.psd .* ...
	  (a^2 + (b + 2*pi*freqLin).^2);
    end
  end
  
  if noiseSourceIsFlicker(noise.contrib(i))
    fprintf(1, 'noise source %s is a 1/f noise source\n', ...
	noise.contrib(i).source.name);
    noise.contrib(i).linout.psd = noise.contrib(i).linout.psd ./ freqLin;
  end
  
  noise.contrib(i).integratedOut = trapz(noise.contrib(i).linout.psd) ;
  noise.contrib(i).averagedOut = noise.contrib(i).integratedOut / freqBand ;
  
  if computeEqIn
    if isfield(noise.tf, 'zeros') & isfield(noise.tf, 'poles')
      noise.contrib(i).lineqIn.psd = noise.contrib(i).linout.psd ./ noise.tf.dcVal^2;
      for j = 1:length(noise.tf.zeros)
    	a = real(noise.tf.zeros(j).val);
    	b = imag(noise.tf.zeros(j).val);
	    noise.contrib(i).lineqIn.psd = (a^2 + b^2) * noise.contrib(i).lineqIn.psd ./ ...
	       (a^2 + (b + 2*pi*freqLin).^2);
      end
      for j = 1:length(noise.tf.poles)
    	a = real(noise.tf.poles(j).val);
    	b = imag(noise.tf.poles(j).val);
    	noise.contrib(i).lineqIn.psd = 1/(a^2 + b^2) * noise.contrib(i).lineqIn.psd .* ...
	     (a^2 + (b + 2*pi*freqLin).^2);
      end
    end
  end

  noise.contrib(i).integratedEqIn = trapz(noise.contrib(i).lineqIn.psd) ;
  noise.contrib(i).averagedEqIn = noise.contrib(i).integratedEqIn / freqBand ;
  
  clear noise.contrib(i).linout.psd ;
  clear noise.contrib(i).lineqIn.psd ;
end

% computation of the total output noise power integrated spectral density:
noise.integratedOut = 0;
noise.averagedOut = 0;
for i = 1:noise.nContribs
  noise.integratedOut = noise.integratedOut + noise.contrib(i).integratedOut;
  noise.averagedOut = noise.averagedOut + noise.contrib(i).averagedOut;
end

% computation of the total input-referred noise power spectral density:
if computeEqIn
  noise.integratedEqIn = 0;
  noise.averagedEqIn = 0;
  for i = 1:noise.nContribs
    noise.integratedEqIn = noise.integratedEqIn + noise.contrib(i).integratedEqIn;
    noise.averagedEqIn = noise.averagedEqIn + noise.contrib(i).averagedEqIn;
  end
end
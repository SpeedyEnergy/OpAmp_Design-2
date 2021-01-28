function result = mosCheckFreq(mos, freqName, operatingFreq)
%MOSCHECKFREQ check whether a given operating frequency is higher than fT
%  (c) IMEC, 2004
%  IMEC confidential 
%

switch freqName
  case {'ft'} % can be extended in the future with 'fmax' if this will be
      % tabulated
    if mos.(freqName) < operatingFreq
      warning('operating frequency %g is larger than %s = %g\n', ...
	  operatingFreq, freqName, mos.(freqName));
      result = 0;
    else
      result = 1;
    end
  otherwise
    error(' specified parameter %s is not a known frequency\n', freqName);
end

function cirWriteParam(fid, simulator, param, value)
%CIRWRITEPARAM write a parameter statement for the given simulator
%
%    cirWriteParam(FID, SIMULATOR, PARAM, VALUE) writes a parameter
%    statement with the appropriate syntax for the given SIMULATOR for the
%    given PARAMETER (should be a string). The VALUE of that parameter
%    should also be given as an argument. The parameter statement is
%    written to the given file identifier FID
%
%    See also fopen, cirWriteParams
%
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

if (strcmpi(simulator, 'hspice') | strcmpi(simulator, 'smartspice'))
  fprintf(fid, '.param %s = %g\n', param, value);
else
  error('The specified simulator is not supported');
end

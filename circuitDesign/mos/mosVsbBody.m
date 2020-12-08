function y = mosVsbBody(mos, vgb, vdb, vov, estimate)
%MOSVSBBODY computation of the source voltage of a MOS transistor with body
%effect
%
%    VSB = mosVsbBody(MOS, VGB, VDB, VOV, ESTIMATE) determines the
%    source-bulk voltage VSB of a MOS transistor by using the given TABLE. 
%    A specified value for the overdrive voltage VOV is given as an argument, 
%    as well as a value for the gate-bulk voltage VGB and drain-bulk voltage VDB. 
%
%    The voltage VSB is determined by fixed-point iteration. A start value
%    ESTIMATE needs to be given for this iteration.
%
%   EXAMPLE:
%    
%      Mn1.vsb = mosVsbBody(Mn1, vdd/2, vdd, 0.2, vdd/2 - 0.2);
%
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%


debug = 0;
if debug
  fprintf(1, 'Calling mosVsbBody with arguments:\n');
  fprintf(1, 'vgb = %g, vdb = %g, vov = %g, estimate = %g\n', vgb, vdb, vov, ...
      estimate)
end

itl = 20; % = maximum number of iterations
relError = 0.005; % specified relative error
nIterations = 0; % counter for the number of iterations 
yprev = estimate;
mos.vgs = vgb - estimate;
mos.vds = vdb - estimate;
mos.vsb = estimate;
y = vgb - vov - mosIntValueWref('vth', mos);
while (abs(y - yprev) > abs(relError*y)) & (nIterations < itl)
yprev = y;
mos.vgs = vgb - yprev;
mos.vds = vdb - yprev;
mos.vsb = yprev;
y = vgb - vov - mosIntValueWref('vth', mos);
nIterations = nIterations + 1;
end
if nIterations == itl
    stringAux = strcat('Max. number of iterations reached in getVsbBody.', ...
        sprintf('\nInitial value for source-bulk voltage is %g\n', estimate), ...
        sprintf('Last value for common source-bulk voltage is %g', y));
    warning(stringAux);
end

  

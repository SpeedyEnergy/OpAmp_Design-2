function y = mosFinVs(mos, vgate, vdrain, vov, estimate)
%MOSFINVS computation of the source voltage of a finFET transistor on FD SOI
%
%    VS = mosFinVs(MOS, VGATE, VDRAIN, VOV, ESTIMATE) determines the
%    source voltage VS (referred to ground) of a finFET transistor. 
%    A specified value for the overdrive voltage VOV is given as an argument, 
%    as well as a value for the gate voltage VGATE and the drain
%    voltage VDRAIN (both referred to ground). 
%
%    The voltage VS is determined by fixed-point iteration. A start value
%    ESTIMATE needs to be given for this iteration.
%
%   EXAMPLE:
%    
%      Mn1.vs = mosFinVs(Mn1, vdd/2, 3/4*vdd, 0.2, vdd/2 - 0.4);
%
%
%
%  (c) IMEC, 2005
%  IMEC confidential 
%


debug = 0;
if debug
  fprintf(1, 'Calling mosFinVs with arguments:\n');
  fprintf(1, 'vgate = %g, vdrain = %g, vov = %g, estimate = %g\n', vgate, vdrain, vov, ...
      estimate)
end

itl = 20; % = maximum number of iterations
relError = 0.005; % specified relative error
nIterations = 0; % counter for the number of iterations 
yprev = estimate;
mos.vgs = vgate - estimate;
tableInCheckRange('vgs', mos.vgs, evalin('base', mos.table));
mos.vds = vdrain - estimate;
tableInCheckRange('vds', mos.vds, evalin('base', mos.table));
mos.vsb = 0;
y = vgate - vov - mosIntValueWref('vth', mos);
while (abs(y - yprev) > abs(relError*y)) & (nIterations < itl)
yprev = y;
mos.vgs = vgate - yprev;
tableInCheckRange('vgs', mos.vgs, evalin('base', mos.table));
mos.vds = vdrain - yprev;
tableInCheckRange('vds', mos.vds, evalin('base', mos.table));
y = vgate - vov - mosIntValueWref('vth', mos);
nIterations = nIterations + 1;
end
if nIterations == itl
    stringAux = strcat('Max. number of iterations reached in getVsbBody.', ...
        sprintf('\nInitial value for source-bulk voltage is %g\n', estimate), ...
        sprintf('Last value for common source-bulk voltage is %g', y));
    warning(stringAux);
end

  

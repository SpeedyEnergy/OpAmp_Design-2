function result = tableValueWref(par, t, lg, vgs, vds, vsb)
%TABLEVALUEWREF retrieving the value of an intrinsic MOS parameter for 
%the reference width, directly from a given table.
%                                                                          
%   RESULT = mosIntValueWref(PARAM, TABLE, LENGTH, VGS, VDS, VSB) returns
%   the value of the MOS operating point parameter PARAM 
%   (specified as a string) for a given TABLE.   
%   PARAM must be exist as a field of the given TABLE. 
%   Possible names for PARAM can be found by running tableDisplay(TABLE).
%
%   EXAMPLE :                                                                
%
%      id = tableValueWref('ids', N, 0.18e-6, 0.7, vdd/2, 0)                            
%  
%   See also tableDisplay, mosIntValue, mosWidth, tableWref
%   
%
% The value of the parameter is computed by multidimensional linear 
% interpolation in the argument TABLE. We programmed the interpolation ourselves.
% This is much more efficient than using the Matlab function "interpn".
%
%  (c) IMEC, 2004
%  IMEC confidential 
%



%if not(mosCheckOpParam(elecParam)) 
%    error('Error in getValueWref: %s is not a known parameter\n', par);
%end


% check whether specified arguments are inside the allowed range:
minVgs = min(t.Input.vgs.init, t.Input.vgs.final);
maxVgs = max(t.Input.vgs.init, t.Input.vgs.final);
minVds = min(t.Input.vds.init, t.Input.vds.final);
maxVds = max(t.Input.vds.init, t.Input.vds.final);
minVsb = min(t.Input.vsb.init, t.Input.vsb.final);
maxVsb = max(t.Input.vsb.init, t.Input.vsb.final);
if vgs < minVgs - eps | vgs > maxVgs + eps | vds < minVds - eps | vds > maxVds ...
      + eps | vsb < minVsb - eps | ...
        vsb > maxVsb + eps | lg < t.Input.lg.init - eps | lg > ...
	t.Input.lg.final + eps
    tableInCheckRange('lg', lg, t);
    tableInCheckRange('vgs', vgs, t);
    tableInCheckRange('vds', vds, t);
    tableInCheckRange('vsb', vsb, t);
    error('out of range');
end

j = min(binSearch(lg, t.Input.lg.array), t.Input.lg.nInputs ...
    - 1); 
k = min(fix((vgs - t.Input.vgs.init)/t.Input.vgs.step) + 1, ...
    t.Input.vgs.nInputs - 1);
m = min(fix((vds - t.Input.vds.init)/t.Input.vds.step) + 1, ...
    t.Input.vds.nInputs - 1);
n = min(fix((vsb - t.Input.vsb.init)/t.Input.vsb.step) + 1, ...
    t.Input.vsb.nInputs - 1);

J = lg - t.Input.lg.array(j);
J1 = t.Input.lg.array(j+1) - lg;
K = vgs - t.Input.vgs.array(k);
K1 = t.Input.vgs.array(k+1) - vgs;
M = vds - t.Input.vds.array(m);
M1 = t.Input.vds.array(m+1) - vds;
N = vsb - t.Input.vsb.array(n);
N1 = t.Input.vsb.array(n+1) - vsb;
%fprintf(1, 'J = %g, j = %d, J1 = %g\n', J, j, J1);
%fprintf(1, 'K = %g, k = %d, K1 = %g\n', K, k, K1);
%fprintf(1, 'M = %g, m = %d, M1 = %g\n', M, m, M1);
%fprintf(1, 'N = %g, n = %d, N1 = %g\n', N, n, N1);

deltaLg = t.Input.lg.array(j + 1) - t.Input.lg.array(j);

result = ((((J1 * t.Table.(par)(j,k,m,n) + J * t.Table.(par)(j+1,k,m,n))* K1 + ...
    (J1*t.Table.(par)(j,k+1,m,n) + J*t.Table.(par)(j+1,k+1,m,n)) * K) * M1 + ...
    ((J1 * t.Table.(par)(j,k,m+1,n) + J * t.Table.(par)(j+1,k,m+1,n)) * K1 + ...
    (J1 * t.Table.(par)(j,k+1,m+1,n) + J * t.Table.(par)(j+1,k+1,m+1,n)) * K) ...
    * M) * N1 ...
    + ...
    (((J1 * t.Table.(par)(j,k,m,n+1) + J * t.Table.(par)(j+1,k,m,n+1))* K1 + ...
    (J1*t.Table.(par)(j,k+1,m,n+1) + J*t.Table.(par)(j+1,k+1,m,n+1)) * K) * M1 + ...
    ((J1 * t.Table.(par)(j,k,m+1,n+1) + J * t.Table.(par)(j+1,k,m+1,n+1)) * K1 + ...
    (J1 * t.Table.(par)(j,k+1,m+1,n+1) + J * t.Table.(par)(j+1,k+1,m+1,n+1)) * K) ...
    * M) * N) / (deltaLg * t.Input.vgs.step * t.Input.vds.step * t.Input.vsb.step);

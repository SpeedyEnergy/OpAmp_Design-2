function vth = mosDiodeVth(lg, vov, table)
% computation of vth of a diode connected MOS transistor without body
% effect.
% This is done with a fixed point iteration
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

debug = 0;
vth = vov + table.Model.vto; % start value;
for i = 1:4
  vth = tableValueWref('vth', table, lg, vov + vth, vov + vth, 0);
  if debug
    fprintf(1, 'iteration #%d, vth = %g\n', i, vth);
  end
end

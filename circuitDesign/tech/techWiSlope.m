function n = techWiSlope(table, lg, vds, vsb)
%TECHWISLOPE slope in weak inversion of a MOS transistor for a given
% channel length, VDS and VSB. 
%
%  TECHWISLOPE(TABLE, LG, VDS, VSB) returns the slope in weak inversion for
%  a MOS transistor characterized by the given TABLE. This slope depends on
%  the channel length LG, the source-bulk voltage (not for FD SOI) and the
%  drain-source voltage.
%  This slope is computed as the drain current divided by the 
%  transconductance, both evaluated in weak inversion (in fact at the
%  minimum inversion level given in the TABLE), and multiplied with the
%  thermal voltage.
%
%  (c) IMEC, 2004
%  IMEC confidential 
%
switch tableType(table)
  case 'n'
    mult = 1;
  case 'p'
    mult = -1;
end

vgs = mult * min(abs(tableInArray('vgs', table)));

k = 1.38e-23;
q = 1.602e-19;
ut = 300 * k/q;

n = tableValueWref('ids', table, lg, vgs, vds, vsb) / tableValueWref('gm', ...
    table, lg, vgs, vds, vsb) / ut; 

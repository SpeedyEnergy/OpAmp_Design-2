function w = mosWidth(elecParam, elecParamValue, mos)
%MOSWIDTH returns MOS transistor width for a given operating point parameter 
%   
%   W  = mosWidth(ELECPARAM, VALUE, MOS) returns the width W of a MOS
%   transistor. This returned width W is such that a wanted value VALUE for
%   the MOS operating point parameter ELECPARAM (specified as a string) is
%   realized. ELECPARAM must be included in the TABLE of a transistor of the same
%   type. Possible names for PARAM can be found by running tableDisplay on
%   mosTable(mos). If the returned width is below the critical width a warning 
%   is given. If the computed width is below the minimum width for the given 
%   technology, then an error is given.
%
%   See also mosTable, tableDisplay, tableWcrit, tableWmin, tableWref
%
%   EXAMPLE :
%
%     Mn1.w = mosWidth('ids', 0.001, Mn1);
%
%     In this example, the width of Mn1 is computed as Mn1.w, based on a
%     specification of 1mA for the drain current of Mn1. It is assumed that the
%     terminal voltages and the length of Mn1 are already determined.
%
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%


table = mosTable(mos);

if not(mosCheckOpParam(elecParam)) 
    error('Error in getW: %s is not a known parameter\n', elecParam);
end

if (ismember(elecParam, tableWIndepParamNames(table)))
  error('Error: %s cannot be used to determine transistor width\n', ...
      elecParam)
end


refVal = tableValueWref(elecParam, mosTable(mos), ...
    mos.lg, mos.vgs, mos.vds, mos.vsb);
if (ismember(elecParam, tableWInverseParamNames(table)))
  w = tableWref(table) * refVal / elecParamValue; 
else
  % parameter scales linearly proportional to W
  w = tableWref(table) * elecParamValue / refVal; 
end

mosCheckWidth(w, table, mos.name);



function gmOverIds = techGmOverIdsMax(table, lg, vds, vsb)
%TECHGMOVERIDSMAX maximum value of gm/Ids from a given table
%
%  GMOVERIDS = TECHGMOVERIDSMAX(TABLE, LG, VDS, VSB) returns the maximum
%  value of gm/Ids for the given value of channel length LG, drain-source
%  voltage VDS and source-bulk voltage VSB of a MOS transistor, that can be
%  obtained from the given TABLE.
%  This is nothing else but gm/Ids evaluated at the minimum VGS value for
%  an n-MOS transistor (at the maximum VGS value for a p-MOS transistor).
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

gmOverIds = tableValueWref('gm', table, lg, vgs, vds, vsb) / tableValueWref('ids', ...
    table, lg, vgs, vds, vsb); 

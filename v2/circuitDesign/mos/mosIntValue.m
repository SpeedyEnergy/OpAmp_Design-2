function result = mosIntValue(elecParam, mos)
%MOSINTVALUE gets the value of a intrinsic MOS parameter 
%
%    PARAMVALUE = mosIntValue(ELECPARAM, MOS) returns the
%    value of a MOS operating point parameter ELECPARAM (specified as a string) for
%    transistor MOS. The value is computed
%    in S.I. units.   
%    ELECPARAM must be included in the table of operating point parameter values
%    of a transistor. Possible names for ELECPARAM can be found by running
%    tableDisplay(mosTable(mos))  or a parameter with name PARAM can be
%   checked on its validity using mosCheckOpParam(PARAM).                                           
%
%   See also tableDisplay, tableValueWref, mosWidth, tableWref, tableWcrit
%   
%
%    EXAMPLE :
%
%      ids = mosIntValue('ids', Mn1)
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%


table = mosTable(mos);

debug = 0;
if debug
  fprintf(1, 'elecParam is %s, lg is %g, vgs is %g, vds is %g, vsb is %g\n', ...
      elecParam, mos.lg, mos.vgs, mos.vds, mos.vsb);
end



if (ismember(elecParam, tableWIndepParamNames(table)))
  result = mosIntValueWref(elecParam, mos);
elseif (ismember(elecParam, tableWInverseParamNames(table)))
  result = mosIntValueWref(elecParam, mos) * tableWref(table) / mos.w;
else
  result = mosIntValueWref(elecParam, mos) * mos.w / tableWref(table);
end

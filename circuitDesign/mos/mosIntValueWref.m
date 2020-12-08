function result = mosIntValueWref(elecParam, mos)
%MOSINTVALUEWREF retrieving value of an intrinsic MOS parameter for 
%a given MOS transistor, assuming its width equals the reference width.
%                                                                          
%   RESULT = mosIntValueWref(PARAM, MOS) returns the value of the 
%   MOS operating point parameter PARAM (specified as a string) for a given MOS
%   transistor.   
%   PARAM must be included in the TABLE of intrinisc values of a transistor
%   of the same type of the given transistor. Possible names for PARAM can be
%   found by running tableDisplay(TABLE) or a parameter with name PARAM can be
%   checked on its validity using mosCheckOpParam(PARAM).
%
%   EXAMPLE :                                                                
%
%      id = mosIntValueWref('ids', Mn1)                            
%  
%   See also tableDisplay, mosIntValue, mosWidth, tableWref,
%   tableValueWref, mosCheckOpParam
%   
%
%
% The value of the parameter is computed by linear interpolation in the 
% argument table.
%
%  (c) IMEC, 2004
%  IMEC confidential 
%


table = mosTable(mos);

if not(mosCheckOpParam(elecParam)) 
    error('Error in getValueWref: %s is not a known parameter\n', elecParam);
end




result = tableValueWref(elecParam, table, mos.lg, mos.vgs, mos.vds, mos.vsb); 


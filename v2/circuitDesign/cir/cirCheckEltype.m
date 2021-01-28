function cirCheckEltype(string)
%CIRCHECKELTYPE check whether the specified string is a valid element type.
%
%  cirCheckEltype(STRING) returns nothing if the string corresponds to a valid
%  element type. Otherwise, an error is signaled. Valid element types are
%  'nmos', 'pmos', 'res' (for resistor), 'cap' (for capacitor), 'v' (for
%  voltage source), 'i' (for current source) and subcircuit ('subckt').
%  These names are case sensitive.
%
%  EXAMPLE :
% 
%      checkEltype('nmos')
%         returns nothing, since the element type exists.
%      checkEltype('NMOS')
%          returns an error.
%
% See also cirInit
%
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%



if not(ismember(string, cirSupportedEltypes))
  error('Unknown circuit element type');
end

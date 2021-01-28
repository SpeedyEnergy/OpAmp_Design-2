function number = tableTechNode(table)
%TABLETECHNODE returns the feature size (as a number in meters) of the ...
%    technology for which the table has been made.
%
%  EXAMPLE:
%  for a table corresponding to a 90 nm technology:
%
%  tableTechNode(N)
%
%  returns 9e-8
%
%  (c) IMEC, 2005
%  IMEC confidential 
%

number = table.Info.tech.node;

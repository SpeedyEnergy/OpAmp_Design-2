function result = cirSupportedEltypes
%CIRSUPPORTEDELTYPES lists the element types that are supported in the set of
%  Matlab functions for circuit sizing
%
%  RESULT = CIRSUPPORTEDELTYPES returns a cell array that contains the element
%  types that are supported.
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

result = {'nmos', 'pmos', 'res', 'cap', 'ind', 'v', 'i', 'subckt'};



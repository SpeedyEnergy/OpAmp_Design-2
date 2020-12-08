function result = mosCheckWidth(width, table, mosName)
%MOSCHECKWIDTH comparison of a transistor width to the critical width
%
%   result = mosCheckWidth(MOSNAME, WIDTH, TABLE) returns 0 and gives a warning when
%   the specified WIDTH of the specified MOS transistor with name MOSNAME is below 
%   the critical width of the given TABLE.
%   This critical width can be obtained with tableWcrit(TABLE)
%
%   See also tableWcrit
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%


if (width < tableWcrit(table))
  result = 0;
  warning(' width of MOS transistor %s = %g  is smaller than critical width', ...
        mosName, width);
else
  result = 1;
end

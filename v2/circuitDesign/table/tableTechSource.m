function name = tableTechSource(table)
%TABLETECHSOURCE returns the name of the company, research center, ... from 
%  which originates technology for which the table has been made.
%
%  (c) IMEC, 2005
%  IMEC confidential 
%

name = table.Info.tech.source;

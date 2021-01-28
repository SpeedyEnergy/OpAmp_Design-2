function string = mosName(MOS)
%MOSNAME returns the name of a given MOS transistor as a string
%
%
%    NAME = mosName(MOS) returns a string that contains the name of the
%    given MOS transistor.
%
%    EXAMPLE :
%  
%      nameMn1 = mosName(Mn1);
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%


string = MOS.name;

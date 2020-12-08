function perimeter = mosAcmPerimeter(acm, lsws, lswg)
%MOSACMPERIMETER computation of the perimeter of a source or drain region
%
%   mosAcmPerimeter(ACM, LSWS, LSWG) returns the perimeter of a diffusion region
%   (= source or drain region). This perimeter is computed according to the area 
%   calculation method ACM with the given perimeter of the diffusion region under
%   the gate (argument LSWG) and not under the gate (argument LSWS).
%
%  (c) IMEC, 2004
%  IMEC confidential 
%


switch acm
  case {2, 4}
    perimeter = lsws + lswg;
  case 3
    perimeter = lsws;
  otherwise
    error('Area calculation method %d is not implemented\n', acm);
end


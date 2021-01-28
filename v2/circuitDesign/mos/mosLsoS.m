function lsos = mosLsoS(wco, eactco, dcopss)
%MOSLSOS computation of length of source diffusion (not between two poly
%    stripes)
%
% lsos = mosLsoS(WCO, EACTCO, DCOPSS) computes the length of the source
% diffusion not between two poly stripes as a function of contact width WCO,
% the extension of active over contact EACTCO, and the distance between contact
% and poly stripe in the source area DCOPSS
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%


lsos = wco + eactco + dcopss;

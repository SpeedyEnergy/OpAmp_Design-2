function lsod = mosLsoD(wco, eactco, dcopsd)
%MOSLSOS computation of length of the drain diffusion (not between two poly
%    stripes)
%
% lsod = mosLsoD(WCO, EACTCO, DCOPSD) computes the length of the drain
% diffusion not between two poly stripes as a function of contact width WCO,
% the extension of active over contact EACTCO, and the distance between contact
% and poly stripe in the drain area DCOPSD
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%


lsod = wco + eactco + dcopsd;

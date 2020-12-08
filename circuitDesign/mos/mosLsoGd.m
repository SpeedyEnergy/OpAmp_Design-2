function lsogd = mosLsoGd(wco, dcopsd)
%MOSLSOGD computation of length of the drain diffusion between two poly
%    stripes
%
% lsogd = mosLsoGd(WCO, DCOPSD) computes the length of the drain
% diffusion between two poly stripes as a function of contact width WCO and
% the distance between contact and poly stripe in the drain area DCOPSD
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%


lsogd = wco + 2*dcopsd;

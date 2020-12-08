function lsogs = mosLsoGs(wco, dcopss)
%MOSLSOGS computation of length of the source diffusion between two poly
%    stripes
%
% lsogs = mosLsoGs(WCO, DCOPSS) computes the length of the source
% diffusion between two poly stripes as a function of contact width WCO and
% the distance between contact and poly stripe in the source area DCOPSS
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%


lsogs = wco + 2*dcopss;

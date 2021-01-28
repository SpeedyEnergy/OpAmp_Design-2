function MOS = mosNsquares(MOS)
%MOSNSQUARES computation of number of squares of source and drain 
% diffusion for the calculation of extrinsic source and drain resistances of a
% MOS transistor. 
%
%   MOS = mosNsquares(MOS) computes the number of squares of source and drain 
% diffusion. These quantities are added as the fields "nrd" (# squares of drain
% diffusion) and "nrs" (idem for source diffusion), respectively, if at least
% these fields do not exist yet. 
% To compute these quantities, the transistor must have the following
% fields: the width <MOS>.w, the length of drain/source diffusion areas between
% two poly stripes (fields <MOS>.lsogd and <MOS>.lsogs, respectively).
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%


if ~isfield(MOS, 'nrd')
  MOS.nrd = (MOS.lsogd/2) / MOS.w;
end
if ~isfield(MOS, 'nrs')
  MOS.nrs = (MOS.lsogs/2) / MOS.w;
end

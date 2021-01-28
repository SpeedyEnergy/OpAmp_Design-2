function MOS = mosOpValues(MOS)
%MOSOPVALUES computation of all operating parameters of a MOS transistor
%
%   MOS = mosOpValues(MOS) computes the value of intrinsic and extrinsic
%   operating point parameters (in S.I. units) of a given transistor MOS. The
%   datastructure of the transistor which is specified as an argument is also
%   returned. However, after this function has returned, all fields of the
%   transistor related to intrinsic and extrinsic operating point parameters
%   have been added to the datastructure of the transistor. Existing extrinsic
%   geometry parameters are NOT recomputed. Also, the
%   values of vgb, vgd and vdb are computed (or updated) from vgs, vds and vsb. 
%                                                                          
%   EXAMPLE :                                                                
%
%       Mn1 = mosOpValues(Mn1)                          
%
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%


if MOS.nFingers == 1
  warning('Transistor %s has only one finger\n', MOS.name);
end

table = mosTable(MOS);

% intrinsic parameters:
names = fieldnames(table.Table);
for i = 1:length(names)
  if mosCheckOpParam(char(names(i)))
    MOS.(char(names(i))) = mosIntValue(char(names(i)), MOS); 
  end
end  

% extrinsic parameters:
MOS = mosNsquares(MOS); % = computation of "nrd" and "nrs"
MOS = mosAreaPerimeter(MOS);
MOS.cdbE = mosJuncap(mosTable(MOS), MOS.vds + MOS.vsb, MOS.ad, MOS.pdtot, ...
    MOS.w);
MOS.csbE = mosJuncap(mosTable(MOS), MOS.vsb, MOS.as, MOS.pstot, ...
    MOS.w);
% total junction capacitors are sum of intrinsic and extrinsic capacitors:
MOS.cdb = MOS.cdbI + MOS.cdbE;
MOS.csb = MOS.csbI + MOS.csbE;

% computation of vgb, vdb and vgd:
MOS.vgb = MOS.vgs + MOS.vsb;
MOS.vgd = MOS.vgs - MOS.vds;
MOS.vdb = MOS.vds + MOS.vsb;

% computation of fT:
MOS.ft = 1/2/pi*MOS.gm/(MOS.cgs + MOS.cgd + MOS.cgb);

% computation of vov:
MOS.vov = MOS.vgs - MOS.vth;

function MOS = mosAreaPerimeter(MOS)
%MOSAREAPERIMETER computation of areas and perimeters of source and drain 
% regions of a MOS transistor.
%
%   MOS = mosAreaPerimeter(MOS) computes 
%    - area of source and drain regions (fields "ad" and "as"), 
%    - length of the sidewall of the source or drain
%      area which is not under the gate (fields "lsd" and "lss"), 
%    - length of the sidewall of the source or drain area which is under the
%      gate (fields "lgd" and "lgs"), 
%    - "pstot", which is the sum of "lss" and "lgs",
%    - "pdtot", which is the sum of "lsd" and "lgd",
%    - "ps" which is either "pstot" or "pstot - w", depending on the model that
%      is used.
%    - "pd", which is either "pdtot" or "pdtot - w", depending on the model
%      that is used. 
%    These quantities are added as fields to the given MOS transistor if they
%    do not exist yet, else they are not recomputed.
%    To compute these quantities, the transistor must have the following
%    fields: the width <MOS>.w, length of the drain/source diffusion areas not
%    between two poly stripes (fields <MOS>.lsod and <MOS>.lsos, respectively),
%    the length of drain/source diffusion areas between two poly stripes
%    (fields <MOS>.lsogd and <MOS>.lsogs, respectively),
%    the amount of fingers (field <MOS>.nFingers) and the multiplicity
%    <MOS>.mult.
%
%    EXAMPLES :
%
%    example 1:
%    Mn1 = mosAreaPerimeter(Mn1)
%
%    example 2:
%    Mn1.w = 20, Mn1.LSOs = 1, Mn1.LSOd = 1, Mn1.LSOgs = 1, Mn1.LSOgd = 1,
%    Mn1.nFingers = 1, Mn1.mult = 1 
%    Mn1 = mosAreaPerimeter(Mn1)
%    returns Mn1 with the following fields computed
%    Mn1.ad = 20, Mn1.lsd = 22, Mn1.lgd = 20, Mn1.as = 20, Mn1.lss = 22, 
%    Mn1.lgs = 20, Mn1.pstot = 42, Mn1.pdtot = 42
%   
%    example 3:    
%    Mn1.w = 20, Mn1.LSOs = 1, Mn1.LSOd = 1, Mn1.LSOgs = 1, Mn1.LSOgd = 1,
%    Mn1.nFingers = 4, Mn1.mult = 1
%    Mn1 = mosAreaPerimeter(Mn1)
%    returns Mn1 with the following fields computed
%    Mn1.ad = 10, Mn1.lsd = 4, Mn1.lgd = 20, Mn1.as = 15, Mn1.lss = 16,
%    Mn1.lgs = 20, Mn1.pstot = 36, Mn1.pdtot = 24
%
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

debug = 0;
Wnom = MOS.w;
LSOs = MOS.lsos;
LSOd = MOS.lsod;
LSOgs = MOS.lsogs;
LSOgd = MOS.lsogd;
FOLD = MOS.nFingers;
MULT = MOS.mult;

if debug
    fprintf(1, 'Wnom = %sm, LSOs = %sm, LSOd = %sm, LSOgs = %sm\n', ...
        eng(Wnom), eng(LSOs), eng(LSOd), eng(LSOgs));
    fprintf(1, 'LSOgd = %sm, FOLD = %d, MULT = %d\n', ...
        eng(LSOgd), FOLD, MULT);
end
x=(FOLD+1)/2;
d=x-mod(x,1);      %number of drains
s=x+mod(x,1);     %number of sources
if MOS.geo == 2 % i.e. more drain regions than source regions
  d=x+mod(x,1);      %number of drains
  s=x-mod(x,1);     %number of sources
end

% There can be a difference between LSOd and LSOgd and between LSOs and LSOgs:
deltad=LSOd-LSOgd;
deltas=LSOs-LSOgs;

% There is a difference in transistor layout parameters if the FOLD factor is
% odd or even: 
if FOLD == 1 
  FLD = 0;   % not folded
else   
  FLD = 1;       % folded
end
odd=mod(FOLD,2);

WSO = Wnom/FOLD;

if debug
    fprintf(1, 'FOLD = %d, FLD = %d, odd = %d, WS0 = %g, deltad = %g, deltas = %g\n', ...
        FOLD, FLD, odd, WSO, deltad, deltas);
    fprintf(1, 'LSOgd = %g, LSOgs = %g, d = %d, s = %d, LSOd = %g,  LSOs = %g\n', ...
        LSOgd, LSOgs, d, s, LSOd, LSOs);
end

% parameters of the drain junction
if ~isfield(MOS, 'ad')
  MOS.ad=FLD*(WSO*(deltad*odd+LSOgd*d))+(1-FLD)*WSO*LSOd;
end
if ~isfield(MOS, 'lsd')
  MOS.lsd=FLD*(WSO*(2*d-FOLD) + 2*(deltad*odd+LSOgd*d))+(1-FLD)*(WSO+2*LSOd);
end
if ~isfield(MOS, 'lgd')
  MOS.lgd=FOLD*WSO;
end
if ~isfield(MOS, 'pdtot')
  MOS.pdtot = MOS.lsd + MOS.lgd;
end
if ~isfield(MOS, 'pd')
  MOS.pd = MOS.pdtot;
end


% parameters of the source junction
if ~isfield(MOS, 'as')
  MOS.as = FLD*(WSO*(deltas*(2-odd)+LSOgs*s))+(1-FLD)*WSO*LSOs;
end
if ~isfield(MOS, 'lss')
  MOS.lss = FLD*(WSO*(2*s-FOLD) + ...
      2*(deltas*(2-odd)+LSOgs*s))+(1-FLD)*(WSO+2*LSOs);
end
if ~isfield(MOS, 'lgs')
  MOS.lgs=FOLD*WSO;
end
if ~isfield(MOS, 'pstot')
  MOS.pstot = MOS.lss + MOS.lgs;
end
if ~isfield(MOS, 'ps')
  MOS.ps = MOS.pstot;
end

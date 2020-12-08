function value = mosExtValue(param, mos)
%MOSEXTVALUE computation of parameters related to the extrinsic 
%part of a MOS transistor.
%
%    value = mosExtValue(NAME, MOS) returns the value of a parameter related
%    to the extrinsic part of a MOS transistor. The NAME of the parameter
%    should be specified as a string. Possible names are
%    cdbE, csbE, as, ad, ps, pd, lsd, lgd, lss, lgs
%
%    EXAMPLE :
%
%       Mn1.cdbE = mosExtValue('cdbE', Mn1);
%
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

debug = 1;
[ad, lsd, lgd, as, lss, lgs] = mosAreaPerimeter(mos.w, mos.lsos, mos.lsod, ... 
    mos.lsogs, mos.lsogd, mos.nFingers, mos.mult);
if debug
    fprintf(1, 'ad = %sm, as = %sm, lsd = %sm, lgd = %sm, lss = %sm, lgs = %sm\n', ...
        eng(ad), eng(as), eng(lsd), eng(lgd), eng(lss), eng(lgs));
end 
    
acm = tableModelParam('acm', mosTable(mos));

switch param
  case 'cdbE'
    value = mosJuncap(mosTable(mos), mos.vds + mos.vsb, ad, lsd, lgd);
  case 'csbE'
    value = mosJuncap(mosTable(mos), mos.vsb, as, lss, lgs);
  case 'as'
    value = as;
  case 'ad'
    value = ad;
  case 'lsd'
    value = lsd;
  case 'lss'
    value = lss;
  case 'lgd'
    value = lgd;
  case 'lgs'
    value = lgs;
  case 'pd'
    value = mosAcmPerimeter(acm, lsd, lgd);
  case 'ps'
    value = mosAcmPerimeter(acm, lss, lgs);
  otherwise    
    error(...
	'parameter %s unknown w.r.t. the extrinsic part of a MOS transistor\n',...
	inputname(1));
end

	
    
    

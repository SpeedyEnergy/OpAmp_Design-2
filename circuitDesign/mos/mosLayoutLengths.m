function mos = mosLayoutLengths(mos, table)
%MOSLAYOUTLENGTHS computation of lengths of the source/drain diffusion between
%   and not between two poly stripes.
%
% MOS = mosLayoutLengths(MOS, TABLE) computes for the given MOS transistor the
% lengths of the source/drain diffusion between and not between two poly
% stripes. This computation is done based on layout rule information, that is
% stored in the given TABLE, that corresponds to the given MOS transistor. The
% fields lsos, lsod, lsogs and lsogd of the given MOS transistor are filled
% with 
% - the value of the length of the source diffusion not between two poly 
% stripes (field "lsos"), 
% - the value of the length of the drain diffusion not between two poly
% stripes (field "lsod"),
% - the value of the length of the source diffusion between two poly 
% stripes (field "lsogs"),
% - the value of the length of the drain diffusion between two poly 
% stripes (field "lsogd")
% Finally, the MOS transistor is returned.
% The following layout rule data should be present in the field table.Info:
% - table.Info.wco: width of a contact
% - table.Info.eactco: extension of active over contact
% - table.Info.dcopss: distance between contact and poly stripe in the source
% area 
% - table.Info.dcopsd: distance between contact and poly stripe in the drain
% area 
%
%  (c) IMEC, 2004
%  IMEC confidential 
%



debug = 0;
mos.lsos = mosLsoS(table.Model.wco, table.Model.eactco, table.Model.dcopss);
mos.lsod = mosLsoS(table.Model.wco, table.Model.eactco, table.Model.dcopsd);
mos.lsogs = mosLsoGs(table.Model.wco, table.Model.dcopss);
mos.lsogd = mosLsoGd(table.Model.wco, table.Model.dcopsd);

if debug
  fprintf(1, '%s.lsos = %s, %s.lsod = %s, %s.lsogs = %s, %s.lsogd = %s\n', ...
      mos.name, eng(mos.lsos), mos.name, eng(mos.lsod), mos.name, ...
      eng(mos.lsogs), mos.name, eng(mos.lsogd));
end


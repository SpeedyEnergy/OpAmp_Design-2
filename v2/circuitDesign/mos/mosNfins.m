function MOS = mosNfins(MOS)
%MOSNFINS computation of the number of fins per finger of a finFET transistor.
%
% MOS = mosNfins(MOS) computes the number of fins for a finFET 
% transistor if the field MOS.nFins does not exist. If the field 
% exists, then nothing is changed. If the field does not exist, the 
% the computed number of fins is stored in the field MOS.nFins. 
% This function can only work if for the MOS transistor the fields w
% and nFingers exist, and for the table of the MOS transistor the fields 
% finWidth and finHeight exist. It is assumed that the total 
% width of one fin is given by wtotalFin = wheight * 2 + wfin.
% The new width is stored in MOS.newWidth.
%
% Note: it could be interesting to do some optimization between 
% nFingers and nFins to achieve the best approximation of w. This 
% is not implemented...
%
% Example: Mn1 = mosNfins(Mn1);
%
%  (c) IMEC, 2005
%  IMEC confidential 
%

if ~isfield(MOS, 'w')
  error('width of transistor %s does not exist', MOS.name);
end

if ~isfield(MOS, 'nFingers')
  error('number of fingers of transistor %s does not exist', MOS.name);
end

table = evalin('base', MOS.table);
finHeight = eval(['table.Info.finHeight']);
finWidth = eval(['table.Info.finWidth']);

if ~isfield(MOS, 'nFins')
    fingerLength = MOS.w / MOS.nFingers ;
    finWidth = finWidth + 2*finHeight ;
    nFins = fingerLength / finWidth ;
    MOS.nFins = round(nFins) ;
    MOS.newWidth = MOS.nFins * finWidth * MOS.nFingers ;
    fprintf('transistor %s: design width: %sm, layout width: %sm\n', ...
	MOS.name, eng(MOS.w),eng(MOS.newWidth));
end

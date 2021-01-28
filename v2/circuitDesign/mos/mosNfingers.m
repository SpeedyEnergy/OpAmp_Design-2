function MOS = mosNfingers(MOS)
%MOSNFINGERS computation of the number of fingers for a MOS transistor.
%
%
% MOS = mosNfingers(MOS) computes the number of
% fingers for a MOS transistor if the field MOS.nFingers does not exist. If the
% field exists, then nothing is changed. If the field does not exist, then the
% computed number of fingers is stored in the field MOS.nFingers. 
% This function can only work if the fields MOS.w and MOS.geo exist. The
% function always returns an even or odd  number of fingers, depending on the
% value of MOS.geo. 
% The function tries to find a fingerwidth between a minimum and a maximum
% width. These widths are determined as follows: if the fields
% MOS.minFingerWidth and MOS.maxFingerWidth exist, then these values are used
% for minimum and maximum width, respectively. If MOS.minFingerWidth does not
% exist, then the value of Wcrit from the table MOS.table is taken. If
% MOS.maxFingerWidth does not exist, then the value maxFingerWidth is taken
% from the base. If the latter does not exist, then an error is given. 
%
% Example:
%
% Mn1 = mosNfingers(Mn1);
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%


if ~isfield(MOS, 'w')
  error('width of transistor %s does not exist', MOS.name);
end
if ~isfield(MOS, 'geo')
  error('geo of transistor %s does not exist', MOS.name);
end

table = evalin('base', MOS.table);

geo = MOS.geo;
% . . maximum/minimum allowed width of one finger
if isfield(MOS,'maxFingerWidth')
    maxFingerWidth = MOS.maxFingerWidth;
else
    maxFingerWidth = evalin('base', 'maxFingerWidth');
end
if isfield(MOS,'minFingerWidth')
    minFingerWidth = MOS.minFingerWidth;
else
    minFingerWidth = table.Info.Wcrit;
end
% . . find nfingers
if ~isfield(MOS,'nFingers')
  switch geo
    case 0
      % just one finger
      MOS.nFingers = 1;
    case {1,2}
      % even number of fingers, drain in the middle for geo=1
      %                         source                      2
      MOS.nFingers = ceil((MOS.w/2)/maxFingerWidth) * 2;
    case 3
      % odd number of fingers
      MOS.nFingers = ceil((MOS.w)/maxFingerWidth);
      if floor(MOS.nFingers/2) == MOS.nFingers/2
	MOS.nFingers = MOS.nFingers + 1;
      end
    otherwise
      fprintf('Error: geo=%g not known\n', geo)
  end
else
  % nothing to do
end
% . . check minimum finger width
nf = MOS.nFingers;
wf = MOS.w / nf;
if wf < minFingerWidth
    fprintf('Warning: fingerwidth for %s %sm < %sm  ==>> nFingers set to 1\n', ...
	MOS.name, eng(wf), eng(minFingerWidth));  
    MOS.nFingers = 1;
    nf = MOS.nFingers;
    wf = MOS.w / nf;
end

function techplots(table, varargin)
%TECHPLOTS generates many plots to characterize an n-MOS or p-MOS in a CMOS 
%    technology, specified by a given TABLE.
%    TECHPLOTS(TABLE) generates plots of
%      - drain current and the small-signal parameters
%      - intrinsic gain gm/gds and gm/gmb
%      - gm/Ids
%      - threshold voltage as a function of channel length and VDS
%      - threshold voltage as a function of VSB (not for FD SOI)
%      - saturation voltage as a function of channel length and the
%         overdrive voltage
%      - intrinsic capacitors as a function of VDS and VGS
%
%    EXAMPLE:
%
%    techplots(P)
%
%
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

if nargin == 2
  gridSpec = varargin{1};
else
  gridSpec = 'minor';
end

if nargin == 3
  lineWidth = varargin{2};
else
  lineWidth = 1;
end

if nargin == 4
  colorSpec = varargin{3};
else
  colorSpec = 'color';
end

if nargin == 5
  if strcmp(colorSpec, 'color')
    colors = varargin{4};
  else
    colors = {'k'};
  end
else
  colors = {'r', 'g', 'b', 'k', 'm'};
end
  
figNumber = 1;

lgMin = tableLmin(table);
vdd = tableMaxVdd(table);
switch tableType(table)
  case 'n'
    mult = 1;
  case 'p'
    mult = -1;
end

vto = tableModelParam('vto', table);
%gridSpec = 'minor';
colors = {'r', 'g', 'b', 'k', 'm'};
%linewidth = 1;
techIds(table, lgMin, tableWref(table), min(10, tableInLength('vgs', table)), ...
    mult*max(abs(vto) - 0.3, tableInInit('vgs', table)), mult*min(abs(vto) ...
    + 0.7, abs(tableMaxVdd(table))), figNumber, colorSpec, lineWidth, ...
    gridSpec, colors); 
figNumber = 6;
techCaps(table, 'cgs', lgMin, tableWref(table), min(15, tableInLength('vds', ...
    table)), 0, mult*abs(tableMaxVdd(table)), figNumber, colorSpec, lineWidth, ...
    gridSpec, colors); 
figNumber = 7;
techCaps(table, 'cgd', lgMin, tableWref(table), min(15, tableInLength('vds', ...
    table)), 0, mult*abs(tableMaxVdd(table)), figNumber, colorSpec, lineWidth, ...
    gridSpec, colors); 
figNumber = 8;
techCaps(table, 'cgb', lgMin, tableWref(table), min(15, tableInLength('vds', ...
    table)), 0, mult*abs(tableMaxVdd(table)), figNumber, colorSpec, lineWidth, ...
    gridSpec, colors); 
figNumber = 9;
techVdsat(table, mult*0.5*abs(tableMaxVdd(table)), min(10, tableInLength('lg', ...
    table)), lgMin, tableInFinal('lg', table), figNumber, colorSpec, ...
    lineWidth, gridSpec, colors); 
figNumber = 10;
techVth(table, min(10, tableInLength('lg', table)), lgMin, tableInFinal('lg', ...
    table), figNumber, colorSpec, lineWidth, gridSpec, colors); 
figNumber = 11;
techGmId(table, 0.5*mult*abs(tableMaxVdd(table)), min(10, tableInLength('lg', ...
    table)), lgMin, tableInFinal('lg', table), figNumber, colorSpec, lineWidth, ...
    gridSpec, colors); 
figNumber = 12;
techVthBody(table, 0.5*mult*abs(tableMaxVdd(table)), min(10, ...
    tableInLength('lg', table)), lgMin, tableInFinal('lg', table), figNumber, ...
    colorSpec, lineWidth, gridSpec, colors);


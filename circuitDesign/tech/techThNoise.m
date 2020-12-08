function techThNoise(table, lg, w, nVgs, vgsInit, vgsFinal, figNumber, ...
    colorspec, linewidth, gridSpec, varargin) 
%TECHTHNOISE plots the power spectral density of the thermal noise of the 
%   drain current of a MOS transistor as a function of VDS for
%   different values of VGS. 
%
%   techThNoise(TABLE, LG, W, NVGS, VGSINIT, VGSFINAL, FIGNUMBER, COLORSPEC,
%   LINEWIDTH, GRIDSPEC, VARARGIN) generates a plot in the figure with
%   number FIGNUMBER of the thermal noise. For a long-channel transistor in
%   saturation, this is approximately 4*kT*2/3*gm. The x-axis of this plot is
%   |VDS|. The above parameters are generated with NVGS values of VGS, ranging
%   between VGSINIT and VGSFINAL. The length and width of the
%   channel are to be specified as an argument (LG and W,
%   respectively). Depending on the value of the argument COLORSPEC it can be
%   color or black-and-white. In the first case COLORSPEC must be equal to
%   'color'. If COLORSPEC has a different value, then a black-and-white plot is
%   generated. The width of the lines in the plot can be specified with a
%   number, which is the argument LINEWIDTH. With the argument GRIDSPEC one can
%   specify whether minor gridlines need to be drawn or not. When GRIDSPEC has
%   the value 'minor', then minor gridlines are generated, else they are not
%   drawn. The function has two optional arguments: the first is a cell array
%   with colors, that specify the colors for the different plots (example:
%   {'r', 'g', 'b'} limits the colors to red, blue and green). The format of
%   the different elements of this cell array is the same as for the
%   specification of colors in a plot command. 
%   The second optional argument is the value of VSB. If this argument is not
%   specified, then VSB is set to zero.
%
%
%
%  (c) IMEC, 2005
%  IMEC confidential 
%



warning off MATLAB:divideByZero;
% = to suppress the warning "Divide by zero" which occurs when dividing by gm
% gmb 

tableInCheckRange('lg', lg, table);
tableInCheckRange('vgs', vgsInit, table);
tableInCheckRange('vgs', vgsFinal, table);


vgsArray = linspace(vgsInit, vgsFinal, nVgs);
vdsArray = tableInArray('vds', table);
% we deliberately make the lengths of vgsArray and vdsArray unequal, in order
% not to mislead the plot function (see "help plot"):
if length(vdsArray) == length(vgsArray)
  vdsArray(length(vgsArray)+1) = vdsArray(length(vgsArray));
end

lineMarkers = {'none', 'o', 'x', 'd', '+', 's'};
lineStyles = {'-', ':', '--', '-.'};
if nargin > 10
  colors = varargin{1};
else
  colors = {'r', 'g', 'b', 'c', 'm', 'k'};
end
if nargin > 11
  vsb = varargin{2}
  tableInCheckRange('vsb', vsb, table);
else
  vsb = 0;
end

for j = 1:length(vgsArray)
  for k = 1:length(vdsArray)
    fig.noise(j,k) = w/tableWref(table) * tableValueWref('di2_id', table, lg, ...
	vgsArray(j), vdsArray(k), vsb);
  end
end

printString = 'di^{2} (thermal)';
unit = 'A^{2}/Hz';
plotName = 'thermal noise';
figure(figNumber);
clf;
h = plot(abs(vdsArray), fig.noise);
set(gcf, 'name', plotName);
for j = 1:length(vgsArray)
  legendArray{j} = strcat('V_{GS}=', sprintf('%1.2f', vgsArray(j)));
  switch colorspec
    case 'color'
      colorVgs = colors{1 + mod(j - 1, length(colors))};
      lineMarkerVgs = 'none';
      lineStyleVgs = lineStyles{1 + mod(floor((j - 1)/length(colors)), ...
	  length(lineStyles))};
    otherwise
      colorVgs = 'k';
      lineMarkerVgs = lineMarkers{1 + mod(j - 1, length(lineMarkers))};
      lineStyleVgs = lineStyles{1 + mod(floor((j - 1)/length(lineStyles)), ...
	  length(lineStyles))};
  end
  set(h(j), 'LineStyle', lineStyleVgs, 'Color', ...
      colorVgs, 'Marker', lineMarkerVgs, 'LineWidth', linewidth);
end
grid on;
switch gridSpec
  case 'minor'
    hp = get(h(1), 'parent');
    set(hp, 'YMinorTick', 'on');
    set(hp, 'XMinorTick', 'on');
    set(hp, 'YMinorGrid', 'on');
    set(hp, 'XMinorGrid', 'on');
end
legend(legendArray);
titleString = strcat(printString, ...
    sprintf(' (V_{SB} = %g V, W = %g%sm, L = %g%sm, %s-MOS of %s)', ...
    vsb, w*1e6, '\mu', lg*1e6, '\mu', tableType(table), ...
    tableTechName(table))); 
title(titleString, 'fontSize', 14);
xlabel('|V_{DS} (V)|', 'fontSize', 14);
ylabel(strcat(printString, ' (', unit, ')'), 'fontSize', 14);



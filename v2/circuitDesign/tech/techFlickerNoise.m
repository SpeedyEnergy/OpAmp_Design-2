function techFlickerNoise(table, lg, w, nVds, vdsInit, vdsFinal, ...
    figNumber, colorspec, linewidth, gridSpec, varargin) 
%TECHCAPS plots the spectral noise density at 1 Hz of the flicker noise ...
%    component of the drain current of a MOS transistor. 
%
%    techFlickerNoise(TABLE, CAPPARAM, LG, W, NVDS, VDSMIN, VDSMAX, 
%    FIGNUMBER, COLORSPEC, LINEWIDTH, GRIDSPEC, VARARGIN) plots the spectral
%    noise density at 1 Hz (units: A^2/Hz) of the flicker noise component of
%    the drain current of a MOS transistor as a function of |VGS| with |VDS| as a 
%    parameter. The plot is generated for a MOS transistor with channel length
%    LG. The transistor is characterized by the given TABLE of intrinsic
%    operating point parameters. The number of VDS values is the argument NVDS.
%    In this way, NVDS values of VDS are taken, which range linearly between an
%    initial value VDSINIT and a final value VDSFINAL.
%    The plot is generated in the figure with number FIGNUMBER. Depending on
%    the value of the argument COLORSPEC it can be color or black-and-white. In
%    the first case COLORSPEC must be equal to 'color'. If COLORSPEC has a
%    different value then a black-and-white plot is generated.
%    The width of the lines in the plot can be specified with a number, which
%    is the argument LINEWIDTH.
%    With the argument GRIDSPEC one can specify whether minor gridlines need to
%    be drawn or not. When GRIDSPEC has the value 'minor', then minor gridlines
%    are generated, else they are not drawn.
%    The function has two optional arguments: the first is a cell array with
%    colors, that specify the colors for the different plots (example: {'r',
%    'g', 'b'} limits the colors to red, blue and green). The format of the
%    different elements of this cell array is the same as for the specification
%    of colors in a plot command.
%    The second optional argument is the value of VSB. If this argument is not
%    specified, then VSB is set to zero.
%
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%


debug = 0;

tableInCheckRange('lg', lg, table);
tableInCheckRange('vds', vdsInit, table);
tableInCheckRange('vds', vdsFinal, table);

vdsArray = linspace(vdsInit, vdsFinal, nVds);
vgsArray = tableInArray('vgs', table);
% we deliberately make the lengths of vgsArray and vdsArray unequal, in order
% not to mislead the plot function (see "help plot"):
if length(vdsArray) == length(vgsArray)
  vgsArray(length(vdsArray)+1) = vgsArray(length(vdsArray));
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

for j = 1:length(vdsArray)
  for k = 1:length(vgsArray)
    fig(j,k) = tableWref(table)/w * tableValueWref('di2_fn', table, lg, ...
	vgsArray(k), vdsArray(j), vsb);
  end
end

if debug
  fprintf(1, 'fig(1,1) = %g, fig(1,2) = %g, fig(2,1) = %g\n', fig(1,1), ...
      fig(1,2), fig(2,2));
end

printString = 'di^2 (flicker noise)';
unit = 'A^{2}/Hz';



% plot command:
figure(figNumber);
clf;
h = semilogy(abs(vgsArray), fig);

if debug
  fprintf(1, 'lineStyles:\n');
  disp(lineStyles);
  
end

for j = 1:length(vdsArray)
  legendArray{j} = strcat('V_{DS}=', sprintf('%1.2f', vdsArray(j)));
  switch colorspec
    case 'color'
      colorVds = colors{1 + mod(j - 1, length(colors))};
      lineMarkerVds = 'none';
      lineStyleVds = lineStyles{1 + mod(floor((j - 1)/length(colors)), ...
	  length(lineStyles))} ;
    otherwise
      colorVds = 'k';
      lineMarkerVds = lineMarkers{1 + mod(j - 1, length(lineMarkers))};
      lineStyleVds = lineStyles{1 + mod(floor((j - 1)/length(lineStyles)), ...
	  length(lineStyles))} ;
  end
  set(h(j), 'LineStyle', lineStyleVds, 'Color', ...
      colorVds, 'Marker', lineMarkerVds, 'LineWidth', linewidth);
  end


legend(legendArray);
titleString = strcat(sprintf(...
    'di^{2}, flicker noise component,  V_{SB} = %g V, W = %g%sm, L = %g%sm, ', ...
    vsb, 1e6*w, '\mu', 1e6*lg, ...
    '\mu'), sprintf(' %s-MOS of %s', tableType(table), tableTechName(table)));
title(titleString, 'fontSize', 14);
xlabel('|V_{GS}| (V)', 'fontSize', 14);
ylabel(sprintf(' di^{2} (%s)', unit), 'fontSize', 14); 
grid on;
switch gridSpec
  case 'minor'
    hp = get(h(1), 'parent');
    set(hp, 'YMinorTick', 'on');
    set(hp, 'XMinorTick', 'on');
    set(hp, 'YMinorGrid', 'on');
    set(hp, 'XMinorGrid', 'on');
end
set(gcf, 'name', '1/f noise');

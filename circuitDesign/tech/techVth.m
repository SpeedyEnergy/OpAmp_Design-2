function techVth(table, nLg, lgInit, lgFinal, figNumber, colorspec, ...
    linewidth, gridSpec, varargin)
%TECHVTH plots the threshold voltage of a MOS transistor as a function of 
%    |VDS| for different channel lengths
%
%    techVth(TABLE, NLG, LGINIT, LGFINAL, FIGNUMBER, COLORSPEC, LINEWIDTH,
%    GRIDSPEC, VARARGIN) generates a plot of |VT| in the
%    figure with number FIGNUMBER. |VT| is plotted as a function of |VDS|
%    with the channel length as a parameter for the given table TABLE of
%    intrinsic operating point parameters of a MOS transistor. 
%    The number of channel length values is the argument NLG. In this way, NLG
%    values of the channel length are taken, which range linearly between an
%    initial value LGINIT and a final value LGFINAL.
%    The plot is generated in the figure with number FIGNUMBER. Depending on
%    the value of the argument COLORSPEC it can be color or black-and-white. In
%    the first case COLORSPEC must be equal to 'color'. If COLORSPEC has a
%    different value, then a black-and-white plot is generated.
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
%  (c) IMEC, 2004
%  IMEC confidential 
%


debug = 0;

tableInCheckRange('lg', lgInit, table);
tableInCheckRange('lg', lgFinal, table);

lgArray = linspace(lgInit, lgFinal, nLg);
lineMarkers = {'none', 'o', 'x', 'd', '+', 's'};
lineStyles = {'-', ':', '--', '-.'};
if nargin > 8
  colors = varargin{1};
else
  colors = {'r', 'g', 'b', 'c', 'm', 'k'};
end
if nargin > 9
  vsb = varargin{2}
  tableInCheckRange('vsb', vsb, table);
else
  vsb = 0;
end


% we assume little dependence of VT on vgs. Therefore, we choose vgs in the
% middle of the range of vgs values:
vgs = 0.5*(tableInFinal('vgs', table) - tableInInit('vgs', table));
vdsArray = tableInArray('vds', table);
% we deliberately make the lengths of vdsArray and lgArray unequal, in order
% not to mislead the plot function (see "help plot"):
if length(vdsArray) == length(lgArray)
  vdsArray(length(lgArray)+1) = vdsArray(length(lgArray));
end


for j = 1:length(lgArray)
  for k = 1:length(vdsArray)
    vthArray(j,k) = tableValueWref('vth', table, lgArray(j), vgs, vdsArray(k), ...
	vsb);
  end
end

  
figure(figNumber);
clf;
h = plot(abs(vdsArray), abs(vthArray));


for j = 1:length(lgArray)
  if lgArray(j) < 1e-6
    scaleFactor = 1e9;
    lgUnit = 'nm';
    formatString = '%3.0f';
  else
    scaleFactor = 1e6;
    lgUnit = '\mum';
    formatString = '%5.2f';
  end
  legendArray{j} = strcat('L = ', sprintf(formatString, lgArray(j)*scaleFactor), ...
      lgUnit);
  switch colorspec
    case 'color'
      colorLg = colors{1 + mod(j - 1, length(colors))};
      lineMarkerLg = 'none';
      lineStyleLg = lineStyles{1 + mod(floor((j - 1)/length(colors)), ...
	  length(lineStyles))};
    otherwise
      colorLg = 'k';
      lineMarkerLg = lineMarkers{1 + mod(j - 1, length(lineMarkers))};
      lineStyleLg = lineStyles{1 + mod(floor((j - 1)/length(lineStyles)), ...
	  length(lineStyles))};
      
  end
  set(h(j), 'LineStyle', lineStyleLg, 'Color', ...
	colorLg, 'Marker', lineMarkerLg, 'LineWidth', linewidth);
  end
  
legend(legendArray);  
xlabel('|V_{DS}|  (V)', 'FontSize', 14);
ylabel('|V_{T}|  (V)', 'FontSize', 14);
grid on;
switch gridSpec
  case 'minor'
    hp = get(h(1), 'parent');
    set(hp, 'YMinorTick', 'on');
    set(hp, 'XMinorTick', 'on');
    set(hp, 'YMinorGrid', 'on');
    set(hp, 'XMinorGrid', 'on');
end

set(gcf, 'name', 'vT(vDS)');
title(sprintf(...
    '|V_{T}| as a function of V_{DS} and L (V_{SB} = %g V, %s-MOS of %s)', vsb, ...
    tableType(table), tableTechName(table)), 'FontSize', 14); 


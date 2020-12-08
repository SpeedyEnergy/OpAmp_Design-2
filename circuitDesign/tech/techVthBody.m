function techVthBody(table, vds, nLg, lgInit, lgFinal, figNumber, colorspec, ...
    linewidth, gridSpec, varargin)
%TECHVTHBODY plots the threshold voltage of a MOS transistor as a function of 
%    |VSB| for different channel lengths
%
%    techVthBody(TABLE, VDS, NLG, LGINIT, LGFINAL, FIGNUMBER, COLORSPEC,
%    LINEWIDTH, GRIDSPEC, VARARGIN) generates a plot of |VT| as a function of
%    |VSB| with a given value of VDS, with the channel length as a
%    parameter. The MOS transistor is characterized by the given table TABLE of
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
%    The function has one optional argument: a cell array with
%    colors, that specify the colors for the different plots (example: {'r',
%    'g', 'b'} limits the colors to red, blue and green). The format of the
%    different elements of this cell array is the same as for the specification
%    of colors in a plot command.
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

debug = 0;

tableInCheckRange('lg', lgInit, table);
tableInCheckRange('lg', lgFinal, table);
tableInCheckRange('vds', vds, table);

lgArray = linspace(lgInit, lgFinal, nLg);
if debug
  fprintf(1, 'lgArray is\n');
  disp(lgArray)
end
lineMarkers = {'none', 'o', 'x', 'd', '+', 's'};
lineStyles = {'-', ':', '--', '-.'};
if nargin > 9
  colors = varargin{1};
else
  colors = {'r', 'g', 'b', 'c', 'm', 'k'};
end


% we assume little dependence of VT on VGS. Therefore, we choose VGS in the
% middle of the range of VGS values:
vgs = 0.5*(tableInFinal('vgs', table) - tableInInit('vgs', table));

% VSB is restricted to tableMaxVdd(table) - VDS:
switch tableType(table)
  case 'n'
    mult = 1;
  case 'p'
    mult = -1;
  otherwise
    error('tableType is %s, which is not a valid type', tableType(table));
end
vsbMax = tableMaxVdd(table) - abs(vds);
arrayAux = abs(tableInArray('vsb', table));
maxVsbIndex = max(find(arrayAux < vsbMax));
vsbArray = mult * arrayAux(1:maxVsbIndex);
% length of vsbArray is made different from length of lgArray, otherwise the
% plot command might be confused. If by accident the two arrays have the same
% length, then an extra element is added to vsbArray by repeating the element
% vsbArray(length(lgArray):
if length(lgArray) == length(vsbArray)
  vsbArray(length(lgArray) + 1) = vsbArray(length(lgArray));
end

if debug
  fprintf(1, 'vsbArray is\n');
  disp(vsbArray);
end

for j = 1:length(lgArray)
  for k = 1:length(vsbArray)
    vthArray(j,k) = tableValueWref('vth', table, lgArray(j), vgs, vds, ...
	vsbArray(k));
  end
end

  
figure(figNumber);
clf;
h = plot(abs(vsbArray), abs(vthArray));


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
xlabel('|V_{SB}|  (V)', 'FontSize', 14);
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

set(gcf, 'name', 'vT(vSB)');
title(sprintf(...
    '|V_{T}| as a function of V_{SB} and L (V_{DS} = %g V, %s-MOS of %s)', vds, ...
    tableType(table), tableTechName(table)), 'FontSize', 14); 


function techGmId(table, vds, nLg, lgInit, lgFinal, figNumber, colorspec, ...
    linewidth, gridSpec, varargin)
%TECHGMID plots gm/IDS of a MOS transistor as a function of the overdrive 
%    voltage VOV. 
%
%    techGmId(TABLE, VDS, NLG, LGINIT, LGFINAL, FIGNUMBER, COLORSPEC,
%    LINEWIDTH, GRIDSPEC, VARARGIN) plots gm/IDS of a MOS transistor as a
%    function of |VOV| = |VGS - VT| and with different channel lengths, for a
%    specified value of VDS. The
%    plot is generated for a MOS transistor characterized by the
%    given TABLE of intrinsic operating point parameters. The number of channel
%    length values is the argument NLG. In this way, NLG values of the channel
%    length are taken, which range linearly between an initial value LGINIT
%    and a final value LGFINAL.
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
%
%  (c) IMEC, 2004
%  IMEC confidential 
%



debug = 0;

tableInCheckRange('lg', lgInit, table);
tableInCheckRange('lg', lgFinal, table);
tableInCheckRange('vds', vds, table);

lgArray = linspace(lgInit, lgFinal, nLg);
arrayAux = tableInArray('vgs', table);
vgsArray = arrayAux(1:length(arrayAux)-1);
% we deliberately make the lengths of vgsArray and lgArray unequal, in order
% not to mislead the plot function (see "help plot"):
if length(lgArray) == length(vgsArray)
  vgsArray(length(lgArray)+1) = vgsArray(length(lgArray));
  vgsArray(length(lgArray)) = 0.5*(vgsArray(length(lgArray)) + vgsArray(length(lgArray)-1));
end

lineMarkers = {'none', 'o', 'x', 'd', '+', 's'};
lineStyles = {'-', ':', '--', '-.'};

if nargin > 9
  colors = varargin{1};
else
  colors = {'r', 'g', 'b', 'c', 'm', 'k'};
end
if nargin > 10
  vsb = varargin{2}
  tableInCheckRange('vsb', vsb, table);
else
  vsb = 0;
end

if strcmpi(tableType(table), 'p')
  nOrP = -1;
else
  nOrP = 1;
end

figure(figNumber);
clf;
hold on;
lgArrayLength = length(lgArray);
dVgs = diff(vgsArray);
vgsArray = nOrP * vgsArray;
for j = 1:lgArrayLength
  for k = 1:length(vgsArray);
    lnId(j,k) = log(tableValueWref('ids', table, lgArray(j), vgsArray(k), vds, vsb));
  end
end

gmOverId = 0*ones(lgArrayLength, length(vgsArray)-1);

for j = 1:lgArrayLength
    gmOverId(j,:) = diff(squeeze(lnId(j,:))) ./ dVgs;
  
  h(j) = plot(vgsArray(2:end), squeeze(gmOverId(j,:)));
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
xlabel('|V_{GS} - V_{T}|  (V)', 'FontSize', 14);
ylabel('g_{m}/I_{DS}  (V^{-1})', 'FontSize', 14);
grid on;
switch gridSpec
  case 'minor'
    hp = get(h(1), 'parent');
    set(hp, 'YMinorTick', 'on');
    set(hp, 'XMinorTick', 'on');
    set(hp, 'YMinorGrid', 'on');
    set(hp, 'XMinorGrid', 'on');
end


set(gcf, 'name', 'gm/IDS');
titleString = strcat('g_{m}/I_{DS}  ', ...
    sprintf('(%s-MOS of %s, V_{SB} = %g V, V_{DS} = %g V)', tableType(table), ...
    tableTechName(table), vsb, vds));
title(titleString, 'fontSize', 14);
hold off;

% plot(VGS(2:end), diff(log(squeeze(N.Table.ids(1,:,20))))./diff(VGS))
% gm/Ids = d ln(Ids) / dVgs


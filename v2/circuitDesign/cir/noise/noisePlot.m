function handle = noisePlot(noise, inout, varargin)
%NOISEPLOT plotting of noise power spectral densities or spot noise figure
%
% HANDLE = noisePlot(NOISE, INOUT) returns a handle or an array of handles to a
% plot of the noise power spectral density at the output or referred to the
% input, or the spot noise figure as a function of frequency. The input
% argument is a structure NOISE, which contains a field "freq" as well as the
% power spectral densities of interest or the spot noise figure. A description
% of the structure NOISE is given in noiseCompute and noiseEnterSource. in, out, inout
% all is default,  0 means only the total, [1 2 3] means that only
% contributions #1, #2 and #3 are plotted. 

debug = 1;

if nargin > 2
  selector = varargin{1};
  if not(isnumeric(selector))
    error('selector of noise contributions should be an array of numbers');
  end
  if length(selector) == 1
    if selector(1) == 0
      contribIndices = [];
    else
      contribIndices = selector;
    end
  else % length(selector) > 1
    contribIndices = sort(selector);
  end
else % no argument was specified
  contribIndices = 0:length(noise.contrib);
end

if nargin > 3
  figNumber = varargin{2};
  figure(figNumber);
end

if nargin > 4
  colorSpec = varargin{3};
else
  colorSpec = 'color';
end

if nargin > 5
  lineWidth = varargin{4};
else
  lineWidth = 1;
end

if nargin > 2
  gridSpec = varargin{5};
  if not(strcmp(gridSpec, 'major')) & not(strcmp(gridSpec, 'minor'))
    error('Illegal value of grid specification');
  end
else
  gridSpec = 'major';
end


switch inout
  case 'eqIn'
    inoutString = 'eqIn';
    titleString = 'Equivalent input noise power spectral density';
  case 'out'
    inoutString = 'out';
    titleString = 'Output noise power spectral density';
  otherwise
    error('argument "inout" is %s, it should be either "eqIn" or "out"', inout);
end

% we will now compose the string that will be evaluated to draw the plot:

if contribIndices(1) == 0
  beginIndex = 2;
  argString = 'noise.freq, noise.(inoutString).psd,';
  legendArray{1} = 'total'; 
else
  beginIndex = 1;
  argString = '';
end

for i = beginIndex:length(contribIndices)
  argString = strcat(argString, ...
      sprintf('noise.freq, noise.contrib(%d).%s.psd,', ...
      contribIndices(i), inoutString));
  legendArray{i} = noise.contrib(contribIndices(i)).source.name;
  
end

argString = argString(1:length(argString)-1);

if debug
  fprintf(1, 'argString is\n%s\n', argString);
end
  
eval(strcat('handle = ', 'loglog(', argString, ');'));
grid on;

if strcmp(gridSpec, 'minor')
   hp = get(handle(1), 'parent');
   set(hp, 'YMinorTick', 'on');
   set(hp, 'XMinorTick', 'on');
   set(hp, 'YMinorGrid', 'on');
   set(hp, 'XMinorGrid', 'on');
end

colors = {'r', 'g', 'b', 'c', 'm', 'k'};
lineMarkers = {'none', 'o', 'x', 'd', '+', 's'};
lineStyles = {'-', ':', '--', '-.'};

for j = 1:length(contribIndices)
  switch colorSpec
    case 'color'
      colorContrib = colors{1 + mod(j - 1, length(colors))};
      lineMarkerContrib = 'none';
      lineStyleContrib = lineStyles{1 + mod(floor((j - 1)/length(colors)), ...
	  length(lineStyles))};
    otherwise
      colorContrib = 'k';
      lineMarkerContrib = lineMarkers{1 + mod(j - 1, length(lineMarkers))};;
      lineStyleContrib = lineStyles{1 + mod(floor((j - ...
	  1)/length(lineMarkers)), length(lineStyles))};
  end
  set(handle(j), 'LineStyle', lineStyleContrib, 'Color', colorContrib, ...
      'Marker', lineMarkerContrib, 'LineWidth', lineWidth);
end

if debug
  for i = 1:length(legendArray)
    fprintf(1, 'legendArray{%d} = %s\n', i, legendArray{i});
  end
end

legend(legendArray);

switch noise.out.type
  case 'v'
    units = 'V^{2}/Hz';
  case 'i'
    units = 'A^{2}/Hz';
end
ylabel(units);
title(titleString);
xlabel('freq. (Hz)');

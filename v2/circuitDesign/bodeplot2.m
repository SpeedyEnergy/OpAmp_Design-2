function [tf, handle] = bodeplot2(varargin)
% bodeplot: plots of magnitude and phase of a transfer function.
% tf1 = bodeplot(fstart, fstop, npoints, dcGain, poles, zeros, fignumber);
%
% example:
%
% tf1 = bodeplot(1, 1e12, 500, loopgain, poles, zeros, 2);
%

fstart = varargin{1};
fstop = varargin{2};
nfreq = varargin{3};
dcVal = varargin{4};
poles = varargin{5};
zeros = varargin{6};
if nargin > 6
  figNumber = varargin{7};
else
  figNumber = 0;
end

f = logspace(log10(fstart), log10(fstop), nfreq);
s = sqrt(-1)*2*pi*f;


tf = dcVal*ones(size(s));
for i = 1:length(poles)
  tf = tf ./ (1 - s./poles(i).val);
end
for i = 1:length(zeros)
  tf = tf .* (1 - s./zeros(i).val);
end

if figNumber
  figure(figNumber);
  clf;
  markerArray = {'+', 'o', '>', '<'};
  colorArray = {'g', 'b', 'k', 'r'};
  pz = [poles zeros];
  legendMag{1} = 'mag';
  legendPhase{1} = 'phase';

  subplot(2,1,1), semilogx(f, 20*log10(abs(tf)));
  xlabel('frequency (Hz)');
  ylabel('gain (dB)');
  grid on;
  hold on
  subplot(2,1,2), semilogx(f, 180/pi*(imag(log(tf))));
  xlabel('frequency (Hz)');
  ylabel('phase (degrees)');
  grid on;
  hold on
  for i = 1:length(pz)
    markerIndex = 1 + mod(fix((i-1)/length(markerArray)), length(markerArray));
    colorIndex = 1 + mod(i-1, length(colorArray));
    pointSpecifier = strcat(markerArray{markerIndex}, colorArray{colorIndex});
    legendMag{i+1} = pz(i).doc;
    legendPhase{i+1} = pz(i).doc;
    fpz(i) = abs(pz(i).val/2/pi);
    tfpzf(i) = interp1(f, tf, fpz(i));
    subplot(2,1,1), semilogx(fpz(i), 20*log10(abs(tfpzf(i))), pointSpecifier);
    subplot(2,1,2), semilogx(fpz(i), 180/pi*imag(log(tfpzf(i))), pointSpecifier);
  end
subplot(2,1,1), legend(legendMag, -1);
subplot(2,1,2), legend(legendPhase, -1);

end

%% Analog Electronics Session 2: cascade amplifier with Miller compensation

%% Adding paths + Loading MOS tables
addpath(genpath('circuitDesign'));
addpath(genpath('functions'));
addpath(genpath('models'));

% clear;
% close all;
% clc;

load ('UMC65_RVT.mat');

%% Initialize everything
designkitName   = 'umc65';
circuitTitle    = 'Analog Design - Session 3';

%Declaration of the circuit components
elementList.nmos = {'Mn1','Mn2'};
elementList.pmos = {'Mp2'};

spec.VDD            = 1.1;
choice.maxFingerWidth = 10e-6;
choice.minFingerWidth = 200e-9;
simulator       ='spectre';
simulFile       = 0;
simulSkelFile   = 0;

analog = cirInit('analog', circuitTitle, 'top', elementList, spec , choice,...
    designkitName, NRVT, PRVT, simulator, simulFile, simulSkelFile);

analog          = cirCheckInChoice(analog, choice);

fprintf('\n--- First Exercise: Miller cap ---\n');
%% EX1: Circuit
disp('                                      ');
disp('       VDD            VDD                 ');
disp('        |              |                  ');
disp('        R1            Mp2                 ');
disp('        |              |                  ');
disp('        |-----+-Cm--Rm-+----+-----+-OUT2  ');
disp('        |     |        |    |     |       ');
disp(' IN1---Mn1    +-------Mn2   CL    RL      ');
disp('        |              |    |     |       ');
disp('        |              |    |     |       ');
disp('       GND            GND  GND   GND      ');

%% EX: Specs
spec.fGBW   = 100e6;        % [Hz] GBW frequency
spec.GBW    = spec.fGBW * 2 * pi; % radians/s
spec.gain   = 80;          % [] voltage gain
spec.gaindB = 20*log10(spec.gain); % 38dB
spec.CL     = 1e-12;       % [F], Load cap
spec.RL     = 5e3;         % [Ohm]

%% Miller Capacitance and resistance
spec.Cm     = 250e-15;    % [F], miller cap
spec.Rm     = 0.0;     % [F], Load cap

%% Translate our specification into MOS parameters
% No Cm and no RM
% % p_inter_node (p1) << p_output_node (p2)
% % GBW = Mn1.gm / Mn2.cgd;
% % p2 = alpa * GBW; PM > 60, alpha = 2.5;
% % p2 = Mn2.gm / spec.CL
% Mn2.gm = 2.5 * spec.GBW * spec.CL;
% Mn1.gm = Mn2.gm/10; %CM < CL

% Cm present and Rm is not
Mn1.gm = (spec.Cm) * spec.GBW;
Mn2.gm = 2.5 * spec.GBW * spec.CL;


%% EX2: Second stage NMOS: Design choices
VOUT        = spec.VDD/2;   % [V], DC output voltage
Mn2.vov     = 0.0;          % [V], overdrive voltage
Mn2.lg      = 400e-9;       % [m], channel length

%% EX2: Second stage NMOS: Implementation
Mn2.vsb = 0;
Mn2.vds = VOUT;
Mn2.vth = tableValueWref('vth', NRVT, Mn2.lg, 0, Mn2.vds, Mn2.vsb);
Mn2.vgs = Mn2.vov + Mn2.vth;
Mn2.w = mosWidth('gm', Mn2.gm, Mn2);
Mn2     = mosNfingers(Mn2);
Mn2     = mosOpValues(Mn2);

%% EX2: Second stage PMOS: Design choices
Mp2.vov     = -0.2;         % [V], overdrive voltage
Mp2.lg      = 400e-9;        % [m], channel length

%% EX2: Second stage PMOS: Implementation;
Mp2.vsb = 0;
Mp2.vds = VOUT - spec.VDD;
Mp2.ids = Mn2.ids + VOUT/spec.RL;
Mp2.vth = tableValueWref('vth', PRVT, Mp2.lg, 0, Mp2.vds, Mp2.vsb);
Mp2.vgs = Mp2.vov + Mp2.vth;
Mp2.w = mosWidth('ids', Mp2.ids, Mp2);
Mp2     = mosNfingers(Mp2);
Mp2     = mosOpValues(Mp2);

%% EX2: First stage NMOS: Design choices
Mn1.vov     = 0.0;         % [V], overdrive voltage
Mn1.lg      = 250e-9;      % [m], channel length

%% EX2: First stage NMOS: Implementation
Mn1.vsb = 0;
Mn1.vds = Mn2.vgs;
Mn1.vth = tableValueWref('vth', NRVT, Mn1.lg, 0, Mn1.vds, Mn1.vsb);
Mn1.vgs = Mn1.vov + Mn1.vth;
Mn1.w   = mosWidth('gm', Mn1.gm, Mn1);
Mn1     = mosNfingers(Mn1);
Mn1     = mosOpValues(Mn1);

%% First Stage Resistance: given stage 1 parameters
R1 = (spec.VDD - Mn1.vds)/Mn1.ids;

%% EX2: Figures of Merit + plot
Av1 = Mn1.gm/(Mn1.gds + 1/R1);
Av2 = Mn2.gm/(Mn2.gds + Mp2.gds + 1/spec.RL);

% pole = 1/(Resistance * Capacitances) = (admitances)/(capacitances)
p1 = (Mn1.gds + 1/R1)/...
    (Mn2.cgs + Mn1.cdd + (spec.Cm + Mn2.cgd)*Av2);

p2 = (Mp2.gds + Mn2.gds +Mn2.gm + 1/spec.RL )/ ...
    (spec.CL +Mp2.cdd+Mn2.cdd);

z1 = 1/((1/Mn2.gm - spec.Rm)*(spec.Cm + Mn2.cgd));

%% what should be CM
% GBW = Mn1.gm / spec.Cm
Cmcalc = Mn1.gm/ spec.GBW;

%% what should be Rm
% Get rid of the zero -> move the zero to infinte
rm1 = 1/Mn2.gm;

% Place the zero on top of the 1st non dominant pole
% z1 = p2
rm2 = 1/(p2*spec.Cm) + 1/Mn2.gm; 

% Recalculating some parameters with the values from LTSpice
p2real = 150e6*2*pi;
gm2real = 3.75e-3;
rm1real = 1/gm2real;
rm2real = 1/(p2real*(spec.Cm+Mn2.cgd)) + 1/gm2real; 

%% Housekeeping code
gainn = Av1*Av2;
if abs(p1)<abs(p2)
    fGBWn = gainn*p1/2/pi;
else 
    fGBWn = gainn*p2/2/pi;
end

fprintf('\n=== Results EX2 ===\n');
fprintf('\nVINT = %gV\n',Mn2.vgs);
fprintf('\nFirst stage: Gain = %g\n',Av1);
fprintf('Second stage: Gain = %g\n1/gdsMn2 = %gOhm, 1/gdsMp2 = %gOhm\n',Av2,1/Mn2.gds,1/Mp2.gds);
fprintf('Total current consumption: %6.2fmA\n',(Mn1.ids+Mp2.ids)/1e-3);

fprintf('\n\t Spec \t\t Actual\n')
fprintf('fGBW: \t %d MHz \t %d MHz\n',...
    spec.fGBW/1e6,round(fGBWn/1e6));
fprintf('gain: \t %g \t\t %g\n',...
    spec.gain,gainn);

freq    = logspace(1,12,1e3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% if control toolbox in Matlab is available
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s   = tf('s');
TF1 = gainn*(1-s*(1/z1))/((1+s/p1)*(1+s/p2));
figure(1);
bode(TF1,freq); grid on;
h = gcr;
setoptions(h,'FreqUnits','Hz');
title('Frequency response cascaded amplifier');
hold all;
figure(2)
TFfb=feedback(TF1,1);
step(TFfb);
title('Step response cascaded amplifier');
hold all

fprintf('\nTransistors in saturation:\n');
if mosCheckSaturation(Mn1)
    fprintf('Mn1:Success\n')
end
if mosCheckSaturation(Mn2)
    fprintf('Mn2:Success\n')
end
if mosCheckSaturation(Mp2)
    fprintf('Mp2:Success\n')
end

%% Print sizes
analog = cirElementsCheckOut(analog); % Update circuit file with 
% transistor sizes
mosPrintSizesAndOpInfo(1,analog); % Print the sizes of the 
% transistors in the circuit file
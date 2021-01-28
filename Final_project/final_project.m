%% Analog Electronics final project

clc; 
%close all; 
clear;

addpath(genpath('../circuitDesign'));
addpath(genpath('../functions'));
addpath(genpath('../models'));


load ('../models/UMC65_RVT.mat');
%% Initializations
designkitName		= 'umc65';
circuitTitle		= 'Analog Design - Project';
elementList.nmos	= {'Mn3','Mn4','Mn6'};
elementList.pmos	= {'Mp1','Mp2','Mp5','Mp7','Mp8'};

choice.maxFingerWidth = 10e-6;
choice.minFingerWidth = 200e-9;

simulator			= 'spectre';
simulFile			= 0;
simulSkelFile		= 0;
spec				= [];
analog			= cirInit('analog',circuitTitle,'top',elementList,spec,choice,...
						designkitName,NRVT,PRVT,simulator,simulFile,simulSkelFile);
analog			= cirCheckInChoice(analog, choice);

spec.VDD = 1.1;

%% Project: circuit
disp('                                                      ');
disp('  VDD          VDD                    VDD             ');
disp('   |             |                      |              ');
disp('  Mp8-+---------Mp7---------------------Mp5            ');
disp('   |--+          |                      |              ');
disp('   |          +--+--+         node 3->  +-----+---OUT  ');
disp('   |          |     |                   |     |        ');
disp('   |    IN1--Mp1   Mp2--IN2             |     |        ');
disp('   |          |     |                   |     |        ');
disp('   |  node 1->|     |<-node 2           |     Cl       ');
disp('   |          |--+  +------+-Cm---Rm----+     |        ');
disp(' Ibias       Mn3-+-Mn4     |            |     |        ');
disp('   |          |     |      +-----------Mn6    |        ');
disp('   |          |     |                   |     |        ');
disp('  GND        GND   GND                 GND   GND       ');

%% AI: Implement your OpAmp according to your handcalculations

%% NMOS
nmos.plotEnable = 0;
% 1: gm/gds in function of Vov
% 2: gm/gds in function of Vgs
% 3: gm/Ids in function of Vov
nmos.option = 1; 

%% PMOS
pmos.plotEnable = 1;
% 1: gm/gds in function of Vov
% 2: gm/gds in function of Vgs
% 3: gm/Ids in function of Vov
pmos.option = 1; 

%% NMOS

%% Get gm/gds in function of Vov
if nmos.plotEnable == 1
    VDD = spec.VDD;
    VDS = 0.01:0.01:VDD;
    VGS = 0.01:0.01:spec.VDD;
    L = [65 80 100 150 200 500 1000]*1e-9;
    Mn6.vds = spec.VDD/2;
    Mn6.vsb = 0;
    result = [];
if nmos.option == 1
    figure
    for ii=1:length(L)
        Mn6.lg = L(ii);
        Mn6.w = 10*Mn6.lg;
        Mn6.vth = tableValueWref('vth',NRVT,Mn6.lg,0,Mn6.vds,Mn6.vsb);
        for jj=VGS%0:0.01:(VGS(end)-Mn6.vth)
            Mn6.vov = jj - Mn6.vth;
            Mn6.vgs = Mn6.vth + Mn6.vov;
            Mn6	= mosNfingers(Mn6);
            Mn6	= mosOpValues(Mn6);
            result = [result Mn6.gm/Mn6.gds];
        end
        hold on
        plot((VGS-Mn6.vth),result);
        result = [];
    end
    xlabel('V_{ov}');
    ylabel('g_{m}/g_{ds}');
    legend('65nm','80nm','100nm','150nm','200nm','500nm','1um');
    title('g_{m}/g_{ds} in function of V_{ov}(w=10L,Vds=VDD/2=0.55V). NMOS');
end

%% Get gm/gds in function of Vgs
if nmos.option == 2
    figure
    for ii=1:length(L)
        Mn6.lg = L(ii);
        Mn6.w = 10*Mn6.lg;
        Mn6.vth = tableValueWref('vth',NRVT,Mn6.lg,0,Mn6.vds,Mn6.vsb);
        for jj=0:0.01:VGS(end)
            Mn6.vgs = jj;
            Mn6.vov = Mn6.vgs - Mn6.vth;
            Mn6	= mosNfingers(Mn6);
            Mn6	= mosOpValues(Mn6);
            result = [result Mn6.gm/Mn6.gds];
        end
        hold on
        plot(0:0.01:VGS(end),result);
        result = [];
    end
    xlabel('V_{GS}');
    ylabel('g_{m}/g_{ds}');
    legend('65nm','80nm','100nm','150nm','200nm','500nm','1um');
    title('g_{m}/g_{ds} in function of V_{GS}(w=10L,Vds=VDD/2=0.55V). NMOS');
end

%% Get gm/Ids in function of Vov
if nmos.option == 3
     figure
    for ii=1:length(L)
        Mn6.lg = L(ii);
        Mn6.w = 10*Mn6.lg;
        Mn6.vth = tableValueWref('vth',NRVT,Mn6.lg,0,Mn6.vds,Mn6.vsb);
        for jj=VGS
            Mn6.vov = jj - Mn6.vth;
            Mn6.vgs = Mn6.vth + Mn6.vov;
            Mn6	= mosNfingers(Mn6);
            Mn6	= mosOpValues(Mn6);
            result = [result Mn6.gm/Mn6.ids];
        end
        hold on
        plot(VGS-Mn6.vth,result);
        result = [];
    end
    xlabel('V_{ov}');
    ylabel('g_{m}/I_{ds}');
    legend('65nm','80nm','100nm','150nm','200nm','500nm','1um');
    title('g_{m}/I_{ds} in function of V_{ov}(w=10L,Vds=VDD/2=0.55V). NMOS');
end
end

%% PMOS

%% Get gm/gds in function of Vov
if pmos.plotEnable == 1
    VDD = spec.VDD;
    VDS = 0.01:0.01:VDD;
    VGS = 0.01:0.01:spec.VDD;
    L = [65 80 100 150 200 500 1000]*1e-9;
    Mp7.vds = -spec.VDD/2;
    Mp7.vsb = 0;
    result = [];
if pmos.option == 1
    figure
    for ii=1:length(L)
        Mp7.lg = L(ii);
        Mp7.w = 10*Mp7.lg;
        Mp7.vth = tableValueWref('vth',PRVT,Mp7.lg,0,Mp7.vds,Mp7.vsb);
        for jj=-VGS%((-1) * (VGS(end) + Mp7.vth)):0.01:0;
            Mp7.vov = jj-Mp7.vth;
            Mp7.vgs = jj;
            Mp7	= mosNfingers(Mp7);
            Mp7	= mosOpValues(Mp7);
            result = [result Mp7.gm/Mp7.gds];
        end
        hold on
        plot(-VGS-Mp7.vth,result);
        result = [];
    end
    xlabel('V_{ov}');
    ylabel('g_{m}/g_{ds}');
    legend('65nm','80nm','100nm','150nm','200nm','500nm','1um');
    title('g_{m}/g_{ds} in function of V_{ov}(w=10L,Vds=-VDD/2=-0.55V). PMOS');
end

%% Get gm/gds in function of Vgs
if pmos.option == 2
    figure
    for ii=1:length(L)
        Mp7.lg = L(ii);
        Mp7.w = 10*Mp7.lg;
        Mp7.vth = tableValueWref('vth',PRVT,Mp7.lg,0,Mp7.vds,Mp7.vsb);
        for jj=((-1) * VGS(end)):0.01:0;
            Mp7.vgs = jj;
            Mp7.vov = Mp7.vgs - Mp7.vth;
            Mp7	= mosNfingers(Mp7);
            Mp7	= mosOpValues(Mp7);
            result = [result Mp7.gm/Mp7.gds];
        end
        hold on
        plot(((-1) * VGS(end)):0.01:0,result);
        result = [];
    end
    xlabel('V_{GS}');
    ylabel('g_{m}/g_{ds}');
    legend('65nm','80nm','100nm','150nm','200nm','500nm','1um');
    title('g_{m}/g_{ds} in function of V_{GS}(w=10L,Vds=-VDD/2=-0.55V). PMOS');
end

%% Get gm/Ids in function of Vov
if pmos.option == 3
     figure
    for ii=1:length(L)
        Mp7.lg = L(ii);
        Mp7.w = 10*Mp7.lg;
        Mp7.vth = tableValueWref('vth',PRVT,Mp7.lg,0,Mp7.vds,Mp7.vsb);
        for jj=((-1) * (VGS(end) + Mp7.vth)):0.01:0;
            Mp7.vov = jj;
            Mp7.vgs = Mp7.vth + Mp7.vov;
            Mp7	= mosNfingers(Mp7);
            Mp7	= mosOpValues(Mp7);
            result = [result Mp7.gm/Mp7.ids];
        end
        hold on
        plot(((-1) * (VGS(end) + Mp7.vth)):0.01:0,result);
        result = [];
    end
    xlabel('V_{ov}');
    ylabel('g_{m}/I_{ds}');
    legend('65nm','80nm','100nm','150nm','200nm','500nm','1um');
    title('g_{m}/I_{ds} in function of V_{ov}(w=10L,Vds=-VDD/2=-0.55V). PMOS');
end
end

fgbw = 28e6;  % Safety factor to account second pole already having an effect
AvdB = 49;
Av = 10^(AvdB/20);
CL = 30e-12;
pd = 2*pi*fgbw/Av;
Mp1.gm = 1.2*2*pi*fgbw*CL/4;

angle = -(-110+180/pi*atan(2*pi*fgbw/pd));
p2 = 1.5*2*pi*fgbw/tan(angle*pi/180);

Mn6.gm = p2*CL;
Mp7.vov = -0.2;
Mn6.lg = 400e-9;%200e-9;
Mp1.vsb = 0;
Mp1.lg = 200e-9;

Mp7.vsb = 0;
Mp8.vsb = 0;
Mp5.vsb = 0;
Mp5.lg = 1e-6;
Mp8.lg = 1e-6;
Mp7.lg = 1e-6;

Mp1.vov = -0.02;%-0.2;
Mn6.vov = 0.1;%0.17;

Vout_wp = (1.1-0.2+0.05)/2;%(1.1-0.2 + 0.17)/2;
Mn6.vsb = 0;
Mn6.vds = Vout_wp;
Mn6.vth = tableValueWref('vth', NRVT, Mn6.lg, 0, Mn6.vds, Mn6.vsb);
Mn6.vgs = Mn6.vov + Mn6.vth;
Mn6.w = mosWidth('gm', Mn6.gm, Mn6);
Mn6 = mosNfingers(Mn6);
Mn6 = mosOpValues(Mn6);

Mp5.vds = Vout_wp - spec.VDD;

Mn4.vsb = 0;
Mn4.vds = Mn6.vgs;
Mn4.vgs = Mn4.vds;
Mn4.lg = 400e-9;%1e-6;
Mn4.vth = tableValueWref('vth', NRVT, Mn4.lg, Mn4.vgs, Mn4.vds, Mn4.vsb);
Mn4.vov = Mn4.vgs-Mn4.vth;

% Estimate Vth with Vds=0, correct later
Mp1.vth = tableValueWref('vth', PRVT, Mp1.lg, 0, 0, Mp1.vsb);
Mp1.vgs = Mp1.vov + Mp1.vth;
Mp1.vgs
Vcm_in_max = Mp1.vgs + Mp7.vov + spec.VDD
Vcm_in_min = Mp1.vgs - Mp1.vov + Mn6.vgs

Vcm = (Vcm_in_max + Vcm_in_min) / 2;

Mp1.vds =  Mn6.vgs - (-Mp1.vgs+Vcm);
Mp1.vth = tableValueWref('vth', PRVT, Mp1.lg, 0, Mp1.vds, Mp1.vsb);
Mp1.vgs = Mp1.vov + Mp1.vth;
Mp1.w = mosWidth('gm', Mp1.gm, Mp1);
Mp1 = mosNfingers(Mp1);
Mp1 = mosOpValues(Mp1);
Mp2 = cirElementCopy(Mp1, Mp2);


Mn4.ids = Mp1.ids;
Mn4.w = mosWidth('ids', Mn4.ids, Mn4);
Mn4 = mosNfingers(Mn4);
Mn4 = mosOpValues(Mn4);
Mn3 = cirElementCopy(Mn4, Mn3);

Mp7.vds = Vcm - Mp1.vgs - spec.VDD;
Mp7.vth = tableValueWref('vth', PRVT, Mp7.lg, 0, Mp7.vds, Mp5.vsb);
Mp7.vgs = Mp7.vov + Mp7.vth;

Mp5.ids = Mn6.ids;
Mp5.vgs = Mp7.vgs;
Mp8.vgs = Mp7.vgs;

Mp5.w = mosWidth('ids', Mp5.ids, Mp5);
Mp5 = mosNfingers(Mp5);
Mp8.w = Mp5.w;
Mp8 = mosNfingers(Mp8);
Mp7.w = Mp8.w*2*Mp1.ids/Mn6.ids;
Mp7 = mosNfingers(Mp7);
Mp5 = mosOpValues(Mp5);

Mp7 = mosOpValues(Mp7);

Mp8.vds = Mp8.vgs;
Mp8 = mosOpValues(Mp8);

% %% AI: Fill out the empty variables required to plot the transfer-function.
% %  meaning of each variable see comment and
% %  location of nodes see line 31 
% 
AvDC1 = Mp1.gm/(Mp1.gds+Mn3.gds);  % DC gain 1st stage
AvDC2 = Mn6.gm/(Mn6.gds+Mp5.gds);  % DC gain 2nd stage
AvDC = AvDC1*AvDC2;
AvdBMade = db(AvDC);
C1    = Mp1.cgd+Mp1.cdb+2*Mn3.cgs+2*Mn3.cgb+Mn3.cgd;  % Capacitance on node 1
G1    = Mn3.gds+Mp1.gds+Mn3.gm;  % Admittance  on node 1
C2    = (1+AvDC2)*(CL/4+Mn6.cgd);  % Capacitance on node 2
G2    = Mn4.gds+Mp2.gds;  % Admittance  on node 2
CLp = CL+Mn6.cdb+Mp5.cdb+Mn4.cgd;
Cn2 = Mn4.cdb+Mn4.cgd+Mp2.cdb+Mp2.cgd;
C3    = CL/4*(CLp+Cn2)+CLp*Cn2;
G3    = Mn6.gm*CL/4;
Vcm_in_min = -Mp1.vgs-Mp1.vdsat+Mn3.vdsat;
Vcm_in_max = spec.VDD+Mp1.vgs+Mp7.vdsat;
%C3    = CL+Mn6.cdb+Mp5.cdb+(1+AvDC2)*CL/4;  % Capacitance on node 3
%G3    = Mp5.gds+Mn6.gds;  % Admittance  on node 3
% 
% %% AI: Set-up Rm, Cc and CL and calculate the zero required for the transfer-fct
% 
spec.Cm = 30e-12/4;
spec.Cl = 30e-12;
spec.Rm = 1/Mn6.gm; 
z1 = inf;
% 
% %% AI: Fill out the empty variables required for the performance summary
Vin_cm_min  = Mp1.vgs - Mp1.vov + Mn6.vgs;   
Vin_cm_max  = Mp1.vgs + Mp7.vov + spec.VDD;
Vout_cm_min = Mn6.vdsat;                         
Vout_cm_max = spec.VDD+Mp8.vdsat;             
Pdiss       = spec.VDD*(Mp8.ids+Mp7.ids+Mp5.ids);
% 
% %% Sanity check (do not modify)
% 
% disp('======================================');
% disp('=      Transistors in saturation     =');
% disp('======================================');
if mosCheckSaturation(Mp1)
	fprintf('\nMp1:Success\n')
end
if mosCheckSaturation(Mp2)
	fprintf('Mp2:Success\n')
end
if mosCheckSaturation(Mn3)
	fprintf('Mn3:Success\n')
end
if mosCheckSaturation(Mn4)
	fprintf('Mn4:Success\n')
end
if mosCheckSaturation(Mp5)
	fprintf('Mp5:Success\n')
end
if mosCheckSaturation(Mn6)
	fprintf('Mn6:Success\n')
end
if mosCheckSaturation(Mp7)
	fprintf('Mp7:Success\n')
end
if mosCheckSaturation(Mp8)
	fprintf('Mp8:Success\n\n')
end


%% Summary of sizes and biasing points (do not modify)

disp('======================================');
disp('=    Sizes and operating points      =');
disp('======================================');
analog = cirElementsCheckOut(analog); % Update circuit file with 
% transistor sizes
mosPrintSizesAndOpInfo(1,analog); % Print the sizes of the 
% transistors in the circuit file
fprintf('IBIAS\t= %6.2fmA\nRm\t= %6.2f Ohm\nCm\t= %6.2fpF\n\n',Mp8.ids/1e-3,spec.Rm,spec.Cm/1e-12);

%% Performance summary (do not modify)

disp('======================================');
disp('=        Performance                 =');
disp('======================================');

fprintf('\nmetrik        \t result\n');
fprintf('Vin,cm,min [mV] \t%.0f\n',Vin_cm_min/1e-3);
fprintf('Vin,cm,max [mV] \t%.0f\n',Vin_cm_max/1e-3);
fprintf('Vout,cm,min [mV] \t%.0f\n',Vout_cm_min/1e-3);
fprintf('Vout,cm,max [mV] \t%.0f\n',Vout_cm_max/1e-3);
fprintf('Pdiss [mW]       \t%.1f\n',Pdiss/1e-3);

%% Ploting transfer function (do not modify)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if control toolbox in Matlab is available
%%%%%%%%%%%%%%s%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s   = tf('s');
% exact transfer function
% TF1 =                     AvDC1*AvDC2*(1+s*C1/(2*G1))*(1-s*(1/z1))/ ...
%      ((1 + s*( spec.Cm*(1+AvDC2)/G2 + (spec.Cm+spec.Cl)/G3 ) + s*s*((spec.Cm*spec.Cl)/(G2*G3)))*(1+s*C1/G1));
TF1 = AvDC1*AvDC2*((1+s*C1/(2*G1))*(1-s*(1/z1)))/ ...
                  ((1+s*C1/G1)*(1+s*C2/G2)*(1+s*C3/G3));


freq = logspace(1,12,1e3);
figure
bode(TF1,2*pi*freq); grid on;
h = gcr;
setoptions(h,'FreqUnits','Hz');
title('Frequency response Opamp');
hold all



%% Analog Electronics final project

clc; 
close all; 
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
nmos.option = 3; 

%% PMOS
pmos.plotEnable = 1;
% 1: gm/gds in function of Vov
% 2: gm/gds in function of Vgs
% 3: gm/Ids in function of Vov
pmos.option = 3; 

%% NMOS

%% Get gm/gds in function of Vov
if nmos.plotEnable == 1
    VDD = spec.VDD;
    VDS = 0.01:0.01:VDD;
    VGS = 0.01:0.01:spec.VDD;
    L = [65 80 100 200 500 1000]*1e-9;
    Mn6.vds = spec.VDD/2;
    Mn6.vsb = 0;
    result = [];
if nmos.option == 1
    figure
    for ii=1:length(L)
        Mn6.lg = L(ii);
        Mn6.w = 10*Mn6.lg;
        Mn6.vth = tableValueWref('vth',NRVT,Mn6.lg,0,Mn6.vds,Mn6.vsb);
        for jj=0:0.01:(VGS(end)-Mn6.vth)
            Mn6.vov = jj;
            Mn6.vgs = Mn6.vth + Mn6.vov;
            Mn6	= mosNfingers(Mn6);
            Mn6	= mosOpValues(Mn6);
            result = [result Mn6.gm/Mn6.gds];
        end
        hold on
        plot(0:0.01:(VGS(end)-Mn6.vth),result);
        result = [];
    end
    xlabel('V_{ov}');
    ylabel('g_{m}/g_{ds}');
    legend('65nm','80nm','100nm','200nm','500nm','1um');
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
    legend('65nm','80nm','100nm','200nm','500nm','1um');
    title('g_{m}/g_{ds} in function of V_{GS}(w=10L,Vds=VDD/2=0.55V). NMOS');
end

%% Get gm/Ids in function of Vov
if nmos.option == 3
     figure
    for ii=1:length(L)
        Mn6.lg = L(ii);
        Mn6.w = 10*Mn6.lg;
        Mn6.vth = tableValueWref('vth',NRVT,Mn6.lg,0,Mn6.vds,Mn6.vsb);
        for jj=0:0.01:(VGS(end)-Mn6.vth)
            Mn6.vov = jj;
            Mn6.vgs = Mn6.vth + Mn6.vov;
            Mn6	= mosNfingers(Mn6);
            Mn6	= mosOpValues(Mn6);
            result = [result Mn6.gm/Mn6.ids];
        end
        hold on
        plot(0:0.01:(VGS(end)-Mn6.vth),result);
        result = [];
    end
    xlabel('V_{ov}');
    ylabel('g_{m}/I_{ds}');
    legend('65nm','80nm','100nm','200nm','500nm','1um');
    title('g_{m}/I_{ds} in function of V_{ov}(w=10L,Vds=VDD/2=0.55V). NMOS');
end
end

%% PMOS

%% Get gm/gds in function of Vov
if pmos.plotEnable == 1
    VDD = spec.VDD;
    VDS = 0.01:0.01:VDD;
    VGS = 0.01:0.01:spec.VDD;
    L = [65 80 100 200 500 1000]*1e-9;
    Mp7.vds = -spec.VDD/2;
    Mp7.vsb = 0;
    result = [];
if pmos.option == 1
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
            result = [result Mp7.gm/Mp7.gds];
        end
        hold on
        plot(((-1) * (VGS(end) + Mp7.vth)):0.01:0,result);
        result = [];
    end
    xlabel('V_{ov}');
    ylabel('g_{m}/g_{ds}');
    legend('65nm','80nm','100nm','200nm','500nm','1um');
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
            Mp7.vov = Mp7.vgs + Mp7.vth;
            Mp7	= mosNfingers(Mp7);
            Mp7	= mosOpValues(Mp7);
            result = [result Mp7.gm/Mp7.gds];
        end
        hold on
        plot(((-1) * VGS(end)):0.01:0,result);
        result = [];
    end
    xlabel('V_{ov}');
    ylabel('g_{m}/g_{ds}');
    legend('65nm','80nm','100nm','200nm','500nm','1um');
    title('g_{m}/g_{ds} in function of V_{ov}(w=10L,Vds=-VDD/2=-0.55V). PMOS');
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
    legend('65nm','80nm','100nm','200nm','500nm','1um');
    title('g_{m}/I_{ds} in function of V_{ov}(w=10L,Vds=-VDD/2=-0.55V). PMOS');
end
end


% %% AI: Fill out the empty variables required to plot the transfer-function.
% %  meaning of each variable see comment and
% %  location of nodes see line 31 
% 
% AvDC1 = ;  % DC gain 1st stage
% AvDC2 = ;  % DC gain 2nd stage
% C1    = ;  % Capacitance on node 1
% G1    = ;  % Admittance  on node 1
% C2    = ;  % Capacitance on node 2
% G2    = ;  % Admittance  on node 2
% C3    = ;  % Capacitance on node 3
% G3    = ;  % Admittance  on node 3
% 
% %% AI: Set-up Rm, Cc and CL and calculate the zero required for the transfer-fct
% 
% spec.Cm = ;
% spec.Cl = ;
% spec.Rm = ; 
% z1 = ;
% 
% %% AI: Fill out the empty variables required for the performance summary
% Vin_cm_min  = ;   
% Vin_cm_max  = ; 
% Vout_cm_min = ;                         
% Vout_cm_max = ;             
% Pdiss       = ;
% 
% %% Sanity check (do not modify)
% 
% disp('======================================');
% disp('=      Transistors in saturation     =');
% disp('======================================');
% if mosCheckSaturation(Mp1)
% 	fprintf('\nMp1:Success\n')
% end
% if mosCheckSaturation(Mp2)
% 	fprintf('Mp2:Success\n')
% end
% if mosCheckSaturation(Mn3)
% 	fprintf('Mn3:Success\n')
% end
% if mosCheckSaturation(Mn4)
% 	fprintf('Mn4:Success\n')
% end
% if mosCheckSaturation(Mp5)
% 	fprintf('Mp5:Success\n')
% end
% if mosCheckSaturation(Mn6)
% 	fprintf('Mn6:Success\n')
% end
% if mosCheckSaturation(Mp7)
% 	fprintf('Mp7:Success\n')
% end
% if mosCheckSaturation(Mp8)
% 	fprintf('Mp8:Success\n\n')
% end
% 
% 
% %% Summary of sizes and biasing points (do not modify)
% 
% disp('======================================');
% disp('=    Sizes and operating points      =');
% disp('======================================');
% analog = cirElementsCheckOut(analog); % Update circuit file with 
% % transistor sizes
% mosPrintSizesAndOpInfo(1,analog); % Print the sizes of the 
% % transistors in the circuit file
% fprintf('IBIAS\t= %6.2fmA\nRc\t= %6.2f Ohm\nCm\t= %6.2fpF\n\n',Mp8.ids/1e-3,spec.Rm,spec.Cm/1e-12);
% 
% %% Performance summary (do not modify)
% 
% disp('======================================');
% disp('=        Performance                 =');
% disp('======================================');
% 
% fprintf('\nmetrik        \t result\n');
% fprintf('Vin,cm,min [mV] \t%.0f\n',Vin_cm_min/1e-3);
% fprintf('Vin,cm,max [mV] \t%.0f\n',Vin_cm_max/1e-3);
% fprintf('Vout,cm,min [mV] \t%.0f\n',Vout_cm_min/1e-3);
% fprintf('Vout,cm,max [mV] \t%.0f\n',Vout_cm_max/1e-3);
% fprintf('Pdiss [mW]       \t%.1f\n',Pdiss/1e-3);
% 
% %% Ploting transfer function (do not modify)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % if control toolbox in Matlab is available
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s   = tf('s');
% % exact transfer function
% TF1 =                     AvDC1*AvDC2*(1+s*C1/(2*G1))*(1-s*(1/z1))/ ...
%      ((1 + s*( spec.Cm*(1+AvDC2)/G2 + (spec.Cm+spec.Cl)/G3 ) + s*s*((spec.Cm*spec.Cl)/(G2*G3)))*(1+s*C1/G1));
% 
% freq = logspace(1,12,1e3);
% figure(1)
% bode(TF1,2*pi*freq); grid on;
% h = gcr;
% setoptions(h,'FreqUnits','Hz');
% title('Frequency response Opamp');
% hold all
% 
% 
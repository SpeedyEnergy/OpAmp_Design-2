% preprocessing of the tables REG018N and REG018P of UMC 0.18um CMOS
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

load('REG018N.mat');
load('REG018P.mat');
global N P;
N = REG018N;
P = REG018P;
clear REG018N REG018P;

techName = 'UMC CMOS 180 nm';
vddMax = 1.8;
modelName = 'Bsim3v3';
wmin = 0.25e-6;
lmin = 0.18e-6;
wcrit = 2e-6;
wref = 1e-5;

modelN.vto = 0.5;
modelN.tox = 4.2e-9;
modelN.cj = 1.03e-3;
modelN.cjsw = 1.34e-10;
modelN.cjswg = modelN.cjsw;
modelN.pb = 0.813;
modelN.pbsw = 0.88;
modelN.pbswg = modelN.pbsw;
modelN.mj = 4.43e-1;
modelN.mjsw = 3.3e-1;
modelN.mjswg = modelN.mjsw;
modelN.hdif = 2.6e-7;
modelN.acm = 3;

N = tablePreprocess(N, wmin, lmin, wref, wcrit, ...
    vddMax, modelName, modelN, techName) ;
clear modelN;

modelP.vto = -0.5;
modelP.tox = 4.2e-9;
modelP.cj = 1.14e-3;
modelP.cjsw = 1.74e-10;
modelP.cjswg = modelP.cjsw;
modelP.pb = 0.762;
modelP.pbsw = modelP.pb;
modelP.pbswg = modelP.pbsw;
modelP.mj = 3.95e-1;
modelP.mjsw = 3.24e-1;
modelP.mjswg = modelP.mjsw;
modelP.hdif = 2.6e-7;
modelP.acm = 3;

P = tablePreprocess(P, wmin, lmin, wref, wcrit, ...
    vddMax, modelName, modelP, techName);
clear modelP;

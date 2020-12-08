load 'UMC18_12';
global N P;
N = N12;
P = P12;
clear N12;
clear P12;
N.Model.mjswg = N.Model.mjsw; 
P.Model.mjswg = P.Model.mjsw; 
N.Model.pbswg = N.Model.pbsw; 
P.Model.pbswg = P.Model.pbsw; 
N.Info.tech.name = 'UMC CMOS 180 nm';
P.Info.tech.name = 'UMC CMOS 180 nm';
N.Info.tech.node = 180e-9;
P.Info.tech.node = 180e-9;
N.Info.tech.source = 'UMC';
P.Info.tech.source = 'UMC';
N.Info.WInverseParams = {'di2_fn'};
P.Info.WInverseParams = {'di2_fn'};
N.Info.Wmin = 0.25e-6;
P.Info.Wmin = 0.25e-6;



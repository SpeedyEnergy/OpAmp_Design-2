fprintf(1, '******\nThis test should run in combination with 180 nm tables\n*******\n');
t0 = clock;
for i = 1:100 
  a.L(i) = tableInterpolNew('ids', N, 0.18e-6+0.1e-6*i, 0.7, 0.949, 0);
end;
interpolTime.L.new = etime(clock,t0);
t0 = clock;

for i = 1:100 
  b.L(i) = tableValueWref('ids', N, 0.18e-6+0.1e-6*i, 0.7, 0.949, 0);
end;
interpolTime.L.old = etime(clock,t0);

interpolError.L = abs((a.L-b.L)./a.L)*100;

t0 = clock;
for i = 1:100 
  a.vgs(i) = tableInterpolNew('ids', N, 0.21e-6, 0.74+0.007*i, 0.949, 0.01);
end;
interpolTime.vgs.new = etime(clock,t0);
t0 = clock;

for i = 1:100 
  b.vgs(i) = tableValueWref('ids', N, 0.21e-6, 0.74+0.007*i, 0.949, 0.01);
end;
interpolTime.vgs.old = etime(clock,t0);

interpolError.vgs = abs((a.vgs-b.vgs)./a.vgs)*100;

t0 = clock;
for i = 1:100 
  a.vsb(i) = tableInterpolNew('ids', N, 0.21e-6, 0.74, 0.949, 0.01+0.002*i);
end;
interpolTime.vsb.new = etime(clock,t0);
t0 = clock;

for i = 1:100 
  b.vsb(i) = tableValueWref('ids', N, 0.21e-6, 0.74, 0.949, 0.01+0.002*i);
end;
interpolTime.vsb.old = etime(clock,t0);

interpolError.vsb = abs((a.vsb-b.vsb)./a.vsb)*100;


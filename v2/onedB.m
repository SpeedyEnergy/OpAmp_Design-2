VoutdB = [-73 -63 -55 -47.56 -40.56 -33.79 -27.11 -20.6 -14.55 -9.72 -7.59 -6.67 -6.34 -6.24 -6.2 -6.19]';
Vin = 10.^(-6:1/3:-1)';

VoutdB = VoutdB+3;
Vin = 2.*Vin;
VindB = db(Vin);
work = [ones(7, 1) VindB(3:9)];
coeffs = work\VoutdB(3:9);
figure
hold on
plot(VindB, VoutdB)
plot(VindB, [ones(length(VindB), 1) VindB]*coeffs)
xlabel('Vin (dBV)')
ylabel('Vout (dB)')
title('Determination of the 1dB-point')

alpha = (VoutdB(10)-VoutdB(9))/(VindB(10)-VindB(9));
p1 = [coeffs(2) coeffs(1)-1];
p2 = [alpha VoutdB(9)-alpha*VindB(9)];
z = roots(p1-p2);
plot([z z], [(p1+[0 1])*[z 1]' p2*[z 1]'], 'gx')
10^(z/20)
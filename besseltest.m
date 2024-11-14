x = logspace(-3,0,500);
Phi = 3 * (sin(x) - x.* cos(x))./x.^3;
q2 = x.^2;
Phi_2 = 1.0 + q2.*(-3./30. + q2.*(3./840. + q2.*(-3./45360.)))
plot(x, Phi-Phi_2,'o')

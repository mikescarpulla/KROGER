T = 800;
P_tot = 1;

x = -0:0.01:5;

n0N2 = 1.5;
n0H2 = 5;
n0NH3 = 0;

GNH3 = G0_NH3_gv(T, P_tot, 1, 'atm');
GN2 = G0_N2_gv(T, P_tot, 1, 'atm');
GH2 = G0_H2_gv(T, P_tot, 1, 'atm');
deltaG = GN2 + 3*GH2 - 2*GNH3
Keq = exp(-deltaG/(8.617e-5*T))

aNH3 = (n0NH3 - 2*x)./((n0NH3 - 2*x)+(n0N2 + x)+(n0H2 + 3*x));
aN2 = (n0N2 + x)./((n0NH3 - 2*x)+(n0N2 + x)+(n0H2 + 3*x));
aH2 = (n0H2 + 3*x)./((n0NH3 - 2*x)+(n0N2 + x)+(n0H2 + 3*x));

Q = aH2.^3 .* aN2 ./ aNH3.^2;


figure(1)
clf
plot(x,Q-Keq)

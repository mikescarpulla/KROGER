mu_O = -8:0.1:8;
P = 1;

figure(1)
clf
hold on
xlim([-10 10])
ylim([-10 10])
xlabel('\mu_{O} (eV)')
ylabel('\mu_{Ga} (eV)')
title('\mu_{Ga} vs. \mu_{O} Constrained by Ga2O3 or GaO, fixed P, different T')


for T = 300:100:2000
[G0_Ga2O3_ZA] = Ga2O3_G0(T,1,'Bar');  % pressure doesnt matter for the solid
[G0_GaO_ZA] = GaO_G0(T,P,'Bar');    %P matters for this as it is a Gas

eq1 = (G0_Ga2O3_ZA - 3*mu_O)/2;
eq2 = G0_GaO_ZA - mu_O;

plot(mu_O,eq1,'r-')
plot(mu_O,eq2,'b-')
input('hit return to plot next T')
end
legend('Ga_2O_3 boundary', 'GaO boundary')



% plot G0 for Ga2O3 and GaO vs T, for different P divide by 10 each time
% starting from 1 atm)
T = 300:100:2000;

figure(2)
clf
hold on
xlim([0 2000])
ylim([-16 2])
xlabel('T (K)')
ylabel('G^o (eV/FU)')

[G0_Ga2O3_ZA] = Ga2O3_G0(T,1,'Bar');
plot(T,G0_Ga2O3_ZA,'r-')

for i=0:-1:-10
    P = 10^i;
    [G0_GaO_ZA] = GaO_G0(T,P,'Bar');

    plot(T,G0_GaO_ZA,'b-')
    input('hit return to plot next P (divide by 10 each time)')
end

legend('G^o_{Ga_2O_3}', 'G^o_{GaO} ')
title('G^o vs T for different GaO partial Pressures ')

%% script to test out the GaN thermodynamic scenarios

conditions.T_equilibrium = 300:10:2000;
conditions.P_tot = 1;
conditions.P_units = 'atm';

conditions.X_GaN = 1;    % this is a solid so Xi=1
conditions.X_Ga = 1;     % set these as the default and change them if needed below

conditions.pN2 = 1;   
conditions.pH2 = 1;
conditions.pNH3 = 1;
%%
figure(1)
clf
hold on

conditions.T_dep_matrix_mu_flag = 'Off';
conditions.mu_GaN = -1.8;   %this is an arbitrary gues, not the real number.  This would be a constant that we get out of a computational paper for the enthalpy of formation for GaN
conditions.mu_conditions_flag = 'Ga-rich';
[mu_Ga, mu_N] = Ga_N2_GaN_chem_potentials_from_conditions(conditions);

plot(conditions.T_equilibrium,mu_Ga,'b-')
plot(conditions.T_equilibrium,mu_N,'g-')

legend('\mu_{Ga} Ga-rich','\mu_{N} Ga-rich')
xlabel('T (K)')
ylabel('\mu (eV)')
ylim([-6 0])
% saveas(gcf, 'Ga-rich.png');
%%
figure(2)
clf
hold on

conditions.mu_conditions_flag = 'N-rich';
conditions.T_dep_matrix_mu_flag = 'Off';
conditions.mu_GaN = -1.8;   %this is an arbitrary gues, not the real number.  This would be a constant that we get out of a computational paper for the enthalpy of formation for GaN
[mu_Ga, mu_N] = Ga_N2_GaN_chem_potentials_from_conditions(conditions);

plot(conditions.T_equilibrium,mu_Ga,'b-')
plot(conditions.T_equilibrium,mu_N,'g-')

legend('\mu_{Ga} N-rich','\mu_{N} N-rich')
xlabel('T (K)')
ylabel('\mu (eV)')
ylim([-6 0])
% saveas(gcf, 'N-rich.png');
%%
figure(3)
clf
hold on

conditions.T_dep_matrix_mu_flag = 'On';
conditions.mu_conditions_flag = 'GaN_touching_Ga(l,s)';

[mu_Ga, mu_N] = Ga_N2_GaN_chem_potentials_from_conditions(conditions);

plot(conditions.T_equilibrium,mu_Ga,'b-')
plot(conditions.T_equilibrium,mu_N,'g-')

legend('\mu_{Ga} GaN touching Ga(l,s)','\mu_{N} GaN touching Ga(l,s)','Location','south')
xlabel('T (K)')
ylabel('\mu (eV)')
ylim([-6 0])
%%
figure(4)
clf
hold on

conditions.mu_conditions_flag = 'Ga_equilibrium_vapor_pressure_over_Ga_at_T';

[mu_Ga, mu_N] = Ga_N2_GaN_chem_potentials_from_conditions(conditions);

plot(conditions.T_equilibrium,mu_Ga,'b-')
plot(conditions.T_equilibrium,mu_N,'g-')

legend('\mu_{Ga} Ga_equilibrium_vapor_pressure_over_Ga_at_T','\mu_{N} Ga_equilibrium_vapor_pressure_over_Ga_at_T')
xlabel('T (K)')
ylabel('\mu (eV)')
ylim([-6 0])
%%
figure(5)
clf
hold on


conditions.T_equilibrium = 300:10:2000;

conditions.mu_conditions_flag = 'pN2-variable';
conditions.P_tot = 1;
conditions.pN2 = 0.8;

conditions.pN2 = conditions.pN2/conditions.P_tot;   %this corresponds to air at 1 atm, or N2+Argon at 1 atm (note that our definition of thermodynamics of Ga-N system does not include reactive oxygen as of now  
% 
% %If we want to do pure N2 at 0.1 atm 
% conditions.P_tot=0.1;
% conditions.X_N2 = 1   %this corresponds to air at 1 atm.


[mu_Ga, mu_N] = Ga_N2_GaN_chem_potentials_from_conditions(conditions);

plot(conditions.T_equilibrium,mu_Ga,'b-')
plot(conditions.T_equilibrium,mu_N,'g-')

legend('\mu_{Ga} pN2-variable','\mu_{N} pN2-variable','Location','south')
xlabel('T (K)')
ylabel('\mu (eV)')
ylim([-6 0])

conditions.P_tot=1;   %
conditions.pN2 = 1;   %set these back to defaults beofre going on
%%
figure(6)
clf
hold on


conditions.mu_conditions_flag = 'N2-1atm';   %this condition definition will include setting X_N2 and Ptot both to 1.  


[mu_Ga, mu_N] = Ga_N2_GaN_chem_potentials_from_conditions(conditions);

plot(conditions.T_equilibrium,mu_Ga,'b-')
plot(conditions.T_equilibrium,mu_N,'g-')

legend('\mu_{Ga} N2-1atm','\mu_{N} N2-1atm','Location','south')
xlabel('T (K)')
ylabel('\mu (eV)')
ylim([-6 0])
%%
figure(7)
clf
hold on

conditions.mu_conditions_flag = 'air-1atm';  

[mu_Ga, mu_N] = Ga_N2_GaN_chem_potentials_from_conditions(conditions);

plot(conditions.T_equilibrium,mu_Ga,'b-')
plot(conditions.T_equilibrium,mu_N,'g-')

legend('\mu_{Ga} air-1atm','\mu_{N} air-1atm','Location','south')
xlabel('T (K)')
ylabel('\mu (eV)')
ylim([-6 0])
%%
figure(8)
clf
hold on



conditions.mu_conditions_flag = 'NH3-H2-N2 Equilibrium';
conditions.P_tot = 0.01;   % this is the final pressure
conditions.n_init_H2 = 0.1;   % these are the INPUT mol fractions of H2 and NH3.  Really it sets the ratio of these two 
conditions.n_init_NH3 = 0.3;
conditions.n_init_N2 = 0.1;
conditions.n_inert = 0.5;

[mu_Ga, mu_N] = Ga_N2_GaN_chem_potentials_from_conditions(conditions);

plot(conditions.T_equilibrium,mu_Ga,'b-')
plot(conditions.T_equilibrium,mu_N,'g-')

legend('\mu_{Ga} N2 determined from NH3/H2 Equilibrium','\mu_{N} N2 determined from NH3/H2 Equilibrium')
xlabel('T (K)')
ylabel('\mu (eV)')
ylim([-6 0])

conditions.P_tot = 1;   % reset to defaults 
conditions.pH2 = 1;   
conditions.pNH3 = 1;
conditions.X_H2 = 1; 
conditions.X_NH3 = 1;
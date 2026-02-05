figure(1)
clf
hold on

P_units = 'atm';
P = 1;


for T = 300:300:2100 %% Plots mu values for 300-2100 K with 300K intervals
    [G0_Ga2O3] = Ga2O3_G0(T,P,P_units);
    
    mu_Ga = G0_Ga2O3/2; %% for O-rich 
    mu_O = G0_Ga2O3/3; %% for Ga-rich 
    

    plot([mu_Ga 0], [0 mu_O])
    xlabel('\mu_{Ga}')
    ylabel('\mu_{O}')
end

legend('300 K', '600 K', '900 K', '1200 K', '1500 K', '1800 K', '2100 K', 'location', 'west')
ax = gca;
ax.YAxisLocation = 'right';
ax.XAxisLocation = 'top';
xlim([-9 0])
ylim([-9 0])


figure(2)
clf
hold on

P_units = 'atm';
P = 1;


for T = 300:300:2100 %% Plots mu values for 300-2100 K with 300K intervals
    [G0_Ga2O] = Ga2O_G0(T,P,P_units);
    
    mu_Ga = G0_Ga2O/2; %% for O-rich 
    mu_O = G0_Ga2O; %% for Ga-rich 
    

    plot([mu_Ga 0], [0 mu_O])
    xlabel('\mu_{Ga}')
    ylabel('\mu_{O}')
end

legend('300 K', '600 K', '900 K', '1200 K', '1500 K', '1800 K', '2100 K', 'location', 'west')
ax = gca;
ax.YAxisLocation = 'right';
ax.XAxisLocation = 'top';
xlim([-9 0])
ylim([-9 0])

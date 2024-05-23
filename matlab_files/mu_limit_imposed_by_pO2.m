function [mu_Si_lim, mu_H_lim, mu_Fe_lim, mu_Sn_lim, mu_Cr_lim, mu_Ti_lim, mu_Ir_lim, mu_Mg_lim, mu_Ca_lim, mu_Zn_lim, mu_Co_lim, mu_Zr_lim, mu_Hf_lim, mu_Ta_lim, mu_Ge_lim] = mu_limit_imposed_by_pO2(T, P_tot, P_units, mu_O)

% order of elements mu_Si mu_H mu_Fe mu_Sn mu_Cr mu_Ti mu_Ir mu_Mg mu_Ca mu_Zn mu_Co mu_Zr mu_Hf mu_Ta mu_Ge


% can make a lext level function on top of this that compares the mu
% requested or computed from constrained element condition to these
% boundaries and tells whether or not a limit is exceeded, and then which
% competing phase is causing the limit.  
% can make a logical vector that says whether the element constraint or the
% competing phase (solubility limit) is limiting the chem potential at each
% T

% for now, starting with the condensed phases only - the G0 for gases
% varies with P, so we need to carefully define the constraints - for
% example that the sum of all partial pressures doesnt exceed the total
% pressure.  But is that actually right?  Let's say the mu of a gas phase
% side product goes huge - that just means it would go away as a gas and
% not stay in the crystal.  


%% all of the condensed phases are called with X_i = 1 for the pure phase.  Only if they alloy with the host phase would we mix the components (e.g. Al2O3 - Ga2O3 )

%% GaOx competing compounds with O %%%%%%%%%%%%%%%%%%%
% Ga2O limit - Ga2O = 2 muGa + muO 
% G0_Ga2O_gv = G0_Ga2O_gv(T, P_tot, 1, P_units);  % gas phase
% Ga_mu_limit_Ga2O = (G0_Ga2O_ls - mu_O)/2;

% GaO limit - GaO = muGa + muO 
% G0_GaO_gv = G0_GaO_gv(T, P_tot, 1, P_units);
% Ga_mu_limit_GaO = G0_Ga2O_ls - mu_O;

% muGa_lim = min([Ga_mu_limit_Ga2O Ga_mu_limit_GaO],[],2);



%% Si oxides   %%%%%%%%%%%%%%%%%%%
% SiO2
G0_SiO2 = G0_SiO2_ls(T, P_tot, 1, P_units);
G0_SiO = G0_SiO_gv(T, P_tot, 1, P_units);
mu_Si_lim = G0_SiO2 - 2*mu_O;

% SiO
% G0_SiO = G0_SiO_ls(T, P_tot, 1, P_units);
% Si_mu_limit_SiO = G0_SiO - mu_O;

% muSi_lim = min([Si_mu_limit_SiO2 Si_mu_limit_SiO],[],2);




%% dihydrogen oxide (water)    %%%%%%%%%%%%%%%%%%%
% G0_H2O_ls = G0_H2O_ls(T, P_tot, 1, P_units);
% G0_H2O_gv = G0_H2O_gv(T, P_tot, 1, P_units);
% muH_lim = (G0_H2O - mu_O)/2;

mu_H_lim = -30*ones(size(T));


%% Fe oxides %%%%%%%%%%%%%%%%%%%
G0_FeO = G0_FeO_ls(T, P_tot, 1, P_units);
G0_Fe2O3 = G0_Fe2O3_ls(T, P_tot, 1, P_units);
G0_Fe3O4 = G0_Fe3O4_ls(T, P_tot, 1, P_units);

Fe_mu_limit_FeO = G0_FeO - mu_O;
Fe_mu_limit_Fe2O3 = (G0_Fe2O3 - 3*mu_O)/2;
Fe_mu_limit_Fe3O4 = (G0_Fe3O4 - 4*mu_O)/3;

mu_Fe_lim = min([Fe_mu_limit_FeO Fe_mu_limit_Fe2O3 Fe_mu_limit_Fe3O4 ],[],2);

figure(8)
clf
hold on
plot(T, Fe_mu_limit_Fe2O3)
plot(T, Fe_mu_limit_Fe3O4)
plot(T, Fe_mu_limit_FeO)
plot(T, mu_Fe_lim,'k--')







%% Sn compounds with O.  List the constraints on mu_Sn from all competing reactions, then take the min allowable chemical potential across all of them.  
% SnO2 limit
G0_SnO2 = G0_SnO2_ls(T, P_tot, 1, P_units);
Sn_mu_limit_SnO2 = G0_SnO2 - 2*mu_O;

% SnO limit
G0_SnO = G0_SnO_ls(T, P_tot, 1, P_units);
Sn_mu_limit_SnO = G0_SnO - mu_O;

mu_Sn_lim = min([Sn_mu_limit_SnO2 Sn_mu_limit_SnO],[],2);







%% Chromium
G0_CrO = G0_CrO_ls(T, P_tot, 1, P_units);
Cr_mu_limit_CrO = G0_CrO - mu_O;

G0_CrO2 = G0_CrO2_ls(T, P_tot, 1, P_units);
Cr_mu_limit_CrO2 = G0_CrO2 - 2*mu_O;

G0_Cr2O3 = G0_Cr2O3_ls(T, P_tot, 1, P_units);
Cr_mu_limit_Cr2O3 = (G0_Cr2O3 - 3*mu_O)/2;

G0_Cr3O4 = G0_Cr3O4_ls(T, P_tot, 1, P_units);
Cr_mu_limit_Cr3O4 = (G0_Cr3O4 - 4*mu_O)/3;

mu_Cr_lim = min([Cr_mu_limit_CrO Cr_mu_limit_CrO2 Cr_mu_limit_Cr2O3 Cr_mu_limit_Cr3O4],[],2);



%% Titanium
G0_TiO = G0_TiO_ls(T, P_tot, 1, P_units);
Ti_mu_limit_TiO = G0_TiO - mu_O;

G0_TiO2 = G0_TiO2_ls(T, P_tot, 1, P_units);
Ti_mu_limit_TiO2 = G0_TiO2 - 2*mu_O;

G0_Ti2O3 = G0_Ti2O3_ls(T, P_tot, 1, P_units);
Ti_mu_limit_Ti2O3 = (G0_Ti2O3 - 3*mu_O)/2;

mu_Ti_lim = min([Ti_mu_limit_TiO Ti_mu_limit_TiO2 Ti_mu_limit_Ti2O3],[],2);




% %% Iridium oxides %%%%%%%%%%%%%%%%%%%
% G0_IrO2_ls = G0_IrO2_ls(T, P_tot, 1, P_units);
% G0_IrO2_gv = G0_IrO2_gv(T, P_tot, 1, P_units);
% G0_IrO3 = G0_IrO3_gv(T, P_tot, 1, P_units);
% 
% muIr_lim_IrO2 = G0_IrO2 - 2*mu_O;
% 
% muIr_lim_IrO2 = G0_IrO2 - 2*mu_O;
% 
% G0_IrO3_gv = G0_IrO2_ls(T, P_tot, 1, P_units);
% 
% 
% 
% 
% muIr_lim = G0_IrO2 - 2*mu_O;
% 

mu_Ir_lim = -30*ones(size(T));



%% MgO   %%%%%%%%%%%%%%%%%%%
G0_MgO = G0_MgO_ls(T, P_tot, 1, P_units);
mu_Mg_lim = G0_MgO - mu_O;



%% CaO  %%%%%%%%%%%%%%%%%%%
G0_CaO = G0_CaO_ls(T, P_tot, 1, P_units);
mu_Ca_lim = G0_CaO - mu_O;


%% ZnO  %%%%%%%%%%%%%%%%%%%
mu_Zn_lim = -30*ones(size(T));    % set to zero if the data is not present - this avoids setting a limit = i.e. you could have pure Zn present


%% Cobalt oxides %%%%%%%%%%%
mu_Co_lim = -30*ones(size(T));    % set to zero if the data is not present - this avoids setting a limit = i.e. you could have pure Zn present


%% ZrO2   %%%%%%%%%%%%%%%%%%%
G0_ZrO2 = G0_ZrO2_ls(T, P_tot, 1, P_units);
mu_Zr_lim = G0_ZrO2 - 2*mu_O;


%% HfO2   %%%%%%%%%%%%%%%%%%%
G0_HfO2 = G0_HfO2_ls(T, P_tot, 1, P_units);
mu_Hf_lim = G0_HfO2 - 2*mu_O;

figure(9)
clf
hold on
plot(T,mu_Hf_lim)


%%% Ta oxides  %%%%%%%%%
mu_Ta_lim = -30*ones(size(T));    % set to zero if the data is not present - this avoids setting a limit = i.e. you could have pure Zn present


%%%%%  GeO, GeO2 %%%% 
mu_Ge_lim = -30*ones(size(T));    % set to -30 if the data is not present - this avoids setting a limit = i.e. you could have pure Zn present



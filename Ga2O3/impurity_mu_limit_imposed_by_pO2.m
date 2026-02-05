function [impurity_mu_lim_out] = impurity_mu_limit_imposed_by_pO2(T, P_tot, P_units, mu_O)
% takes in T, P_tot, and mu_O and figures out the mu values for impurities
% set by their equilibrium with oxygen.  Inherent here is the assumption
% that the impurity concentrations are small, so that they will not form
% compounds with the matrix, for example, or with eachother.  

% To do: improve this to take in any mu conditions (like Ga rich) and find boundaries on
% all the impurities.  

% order of elements Si H Fe Sn Cr Ti Ir Mg Ca Zn Co Zr Hf Ta Ge Pt Rh
% these get transposed then concatenated horizontally into a matrix at end called
% mu_lim_out


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


%% Si oxides   %%%%%%%%%%%%%%%%%%%
% SiO2
G0_SiO2 = G0_SiO2_ls(T, P_tot, 1, P_units);
Si_mu_lim_SiO2 = G0_SiO2 - 2*mu_O;

% SiO
G0_SiO = G0_SiO_gv(T, P_tot, 1, P_units);
Si_mu_lim_SiO = G0_SiO - mu_O;

Si_mu_lim = min([Si_mu_lim_SiO2 Si_mu_lim_SiO],[],2);



%% dihydrogen oxide (water)    %%%%%%%%%%%%%%%%%%%
G0_H2O = min([G0_H2O_ls(T, P_tot, 1, P_units) G0_H2O_gv(T, P_tot, 1, P_units)],[],2);  % take min of G0 2 phases of same substance first.  Then min of mu limits from substaces next. 
H_mu_lim = (G0_H2O - mu_O)/2;



%% Fe oxides %%%%%%%%%%%%%%%%%%%
G0_FeO = G0_FeO_ls(T, P_tot, 1, P_units);
G0_Fe2O3 = G0_Fe2O3_ls(T, P_tot, 1, P_units);
G0_Fe3O4 = G0_Fe3O4_ls(T, P_tot, 1, P_units);

Fe_mu_limit_FeO = G0_FeO - mu_O;
Fe_mu_limit_Fe2O3 = (G0_Fe2O3 - 3*mu_O)/2;
Fe_mu_limit_Fe3O4 = (G0_Fe3O4 - 4*mu_O)/3;

Fe_mu_lim = min([Fe_mu_limit_FeO Fe_mu_limit_Fe2O3 Fe_mu_limit_Fe3O4 ],[],2);  % here an example where mutliple phases could be limiting, find the one setting the most stringent constaint at each T individually 


%% Sn compounds with O.  List the constraints on Sn_mu from all competing reactions, then take the min allowable chemical potential across all of them.  
% SnO2 limit
G0_SnO2 = G0_SnO2_ls(T, P_tot, 1, P_units);
Sn_mu_lim_SnO2 = G0_SnO2 - 2*mu_O;

% SnO limit
G0_SnO = G0_SnO_ls(T, P_tot, 1, P_units);
Sn_mu_lim_SnO = G0_SnO - mu_O;

Sn_mu_lim = min([Sn_mu_lim_SnO2 Sn_mu_lim_SnO],[],2);



%% Chromium
G0_CrO = G0_CrO_ls(T, P_tot, 1, P_units);
Cr_mu_limit_CrO = G0_CrO - mu_O;

G0_CrO2 = G0_CrO2_ls(T, P_tot, 1, P_units);
Cr_mu_limit_CrO2 = G0_CrO2 - 2*mu_O;

G0_Cr2O3 = G0_Cr2O3_ls(T, P_tot, 1, P_units);
Cr_mu_limit_Cr2O3 = (G0_Cr2O3 - 3*mu_O)/2;

G0_Cr3O4 = G0_Cr3O4_ls(T, P_tot, 1, P_units);
Cr_mu_limit_Cr3O4 = (G0_Cr3O4 - 4*mu_O)/3;

Cr_mu_lim = min([Cr_mu_limit_CrO Cr_mu_limit_CrO2 Cr_mu_limit_Cr2O3 Cr_mu_limit_Cr3O4],[],2);



%% Titanium
G0_TiO = G0_TiO_ls(T, P_tot, 1, P_units);
Ti_mu_limit_TiO = G0_TiO - mu_O;

G0_TiO2 = G0_TiO2_ls(T, P_tot, 1, P_units);
Ti_mu_limit_TiO2 = G0_TiO2 - 2*mu_O;

G0_Ti2O3 = G0_Ti2O3_ls(T, P_tot, 1, P_units);
Ti_mu_limit_Ti2O3 = (G0_Ti2O3 - 3*mu_O)/2;

Ti_mu_lim = min([Ti_mu_limit_TiO Ti_mu_limit_TiO2 Ti_mu_limit_Ti2O3],[],2);



% %% Iridium oxides %%%%%%%%%%%%%%%%%%%
G0_IrO2 = min([G0_IrO2_ls(T, P_tot, 1, P_units) G0_IrO2_gv(T, P_tot, 1, P_units)],[],2);
G0_IrO3 = G0_IrO3_gv(T, P_tot, 1, P_units);

Ir_mu_lim_IrO2 = G0_IrO2 - 2*mu_O;
Ir_mu_lim_IrO3 = G0_IrO3 - 3*mu_O;

Ir_mu_lim = min([Ir_mu_lim_IrO2 Ir_mu_lim_IrO3],[],2);



%% MgO   %%%%%%%%%%%%%%%%%%%
G0_MgO = min([G0_MgO_ls(T, P_tot, 1, P_units) G0_MgO_gv(T, P_tot, 1, P_units)],[],2);
Mg_mu_lim = G0_MgO - mu_O;


%% CaO  %%%%%%%%%%%%%%%%%%%
G0_CaO = min([G0_CaO_ls(T, P_tot, 1, P_units) G0_CaO_gv(T, P_tot, 1, P_units)],[],2);
Ca_mu_lim = G0_CaO - mu_O;


%% ZnO  %%%%%%%%%%%%%%%%%%%
Zn_mu_lim = -30*ones(numel(T),1);    % set to zero if the data is not present - this avoids setting a limit = i.e. you could have pure Zn present


%% Cobalt oxides %%%%%%%%%%%
Co_mu_lim = -30*ones(numel(T),1);   


%% ZrO2   %%%%%%%%%%%%%%%%%%%
G0_ZrO2 = min([G0_ZrO2_ls(T, P_tot, 1, P_units) G0_ZrO2_gv(T, P_tot, 1, P_units)],[],2);
Zr_mu_lim = G0_ZrO2 - 2*mu_O;


%% HfO2   %%%%%%%%%%%%%%%%%%%
G0_HfO2 = G0_HfO2_ls(T, P_tot, 1, P_units);
Hf_mu_lim = G0_HfO2 - 2*mu_O;


%%% Ta oxides  %%%%%%%%%
Ta_mu_lim = -30*ones(numel(T),1);    % set to zero if the data is not present - this avoids setting a limit = i.e. you could have pure Zn present


%%%%%  GeO, GeO2 %%%% 
Ge_mu_lim = -30*ones(numel(T),1);    % set to -30 if the data is not present - this avoids setting a limit = i.e. you could have pure Zn present


%%%% Pt oxides  %%%%%%%%
Pt_mu_lim = -30*ones(numel(T),1);    % set to -30 if the data is not present - this avoids setting a limit = i.e. you could have pure Zn present


%%%%% Rh oxides %%%%%%
Rh_mu_lim = -30*ones(numel(T),1);    % set to -30 if the data is not present - this avoids setting a limit = i.e. you could have pure Zn present



%% concatenate them all together
impurity_mu_lim_out = [Si_mu_lim H_mu_lim Fe_mu_lim Sn_mu_lim Cr_mu_lim Ti_mu_lim Ir_mu_lim Mg_mu_lim Ca_mu_lim Zn_mu_lim Co_mu_lim Zr_mu_lim Hf_mu_lim Ta_mu_lim Ge_mu_lim Pt_mu_lim Rh_mu_lim];

if sum(sum(isinf(impurity_mu_lim_out)))
    error('Infinite mu detected - this happens when there is no model for G0 for the T requested for some substance(s).  Figure out which element is causing it and what to do')
end
function [output_holder] = CdTe_fit_defect_energies(defects)


% %% clear workspace and figures   %%%%%%%%%%%%
% clear all
% clc

%% set Physical Constants  %%%%%%%%%
conditions.q = 1.602176565e-19;
conditions.h = 6.62606957e-34;
conditions.kB = 8.6173324e-5;
conditions.mo = 9.1093837e-31;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


T_fit = 1360 : -20 : 900;  % evenly spaced temperatures

% choose upper and lower limits to deltaE to add to each defect's Eform
delta_E_hi = 5;
delta_E_lo = -5;




%% Input the solidus line data - off stoichiometry boundaries for the CdTe solid phase
%% from Greenberg, J Crystal Growth 1996
%% Cd and Te compositions for the boundaries of the solidus region and the mu values for liquidus
%% T in Kelvin.  X from Greenberg originally in %, convert to 0<X<1.  Then to #/cm3 atoms
%% some extra points added to make both Cd and Te data span the whole T range and make them go to 50/50 at low T.

T_Cd_rich = flip([600 697.2 819.6 864.1 871.3 900.8 939.7 969.4 992.0 1073 1123 1173 1223 1300 1360]);    % list of measured points, changed to be decreasing T
pct_Cd_rich_solidus = flip([50 50.0002 50.0001 49.9990 49.9983 49.9990 49.9994 49.9991 49.9967 49.9966 49.9946 49.9957  49.9967 50 50.001]);
% pct_Cd_rich_errs = flip([ 0.0002 0.0001 0.0001 0.0006 0.0001 0.0002 0.0002 0.0003 0.0003 0.0003 0.0003 0.0003]);
% T_Cd_rich = T_Cd_rich(1:end-1);  % delete the lowest-T point as we dont have liquidus data for it
% pct_Cd_rich_solidus = pct_Cd_rich_solidus(1:end-1);
% pct_Cd_rich_errs = pct_Cd_rich_errs(1:end-1);
X_Cd_rich_solidus = pct_Cd_rich_solidus'/100;   % this is now a fraction
% X_Cd_rich_solidus_errs = pct_Cd_rich_errs'/100;    % this is now a fraction

%same for Te
T_Te_rich = flip([600 700 800 900 945.5 1016.6 1046.9 1073 1076.0 1094.1 1103.3 1121.0 1123 1173 1223 1243.1 1282.2 1285.5 1301.0 1312.0 1359.4 1361.8]);
pct_Te_rich_solidus = flip([50 50 50.0002 50.0005000 50.0006 50.0036 50.0055 50.0071 50.0073 50.0092 50.0098 50.0138 50.0138 50.0135 50.0082 50.0087 50.0050 50.0047 50.0044 50.0038 50.0034 50.0013]);
% pct_Te_rich_errs =    flip([ 0.0001  0.0001  0.0004  0.0004  0.0001  0.0003  0.0001  0.0007  0.0007  0.0007  0.0007  0.0008  0.0010  0.0015  0.0030  0.0040  0.0006  0.0004]);
X_Te_rich_solidus = pct_Te_rich_solidus'/100;  % this is now in #/cm3
% X_Te_rich_solidus_errs = pct_Te_rich_errs'/100;  % this is now in #/cm3

%% generate polynomial fits of the data that we can use for fitting calcs to.
X_Cd_poly = polyfit(T_Cd_rich, X_Cd_rich_solidus, 5);
X_Te_poly = polyfit(T_Te_rich, X_Te_rich_solidus, 11);
X_Cd_poly_fit = polyval(X_Cd_poly, T_fit);
X_Te_poly_fit = polyval(X_Te_poly, T_fit);

% plot the modified data and polynomial fits
figure(1)
clf
hold on
plot(T_Cd_rich, X_Cd_rich_solidus, 'bo')  % Greenberg data
plot(T_Te_rich, X_Te_rich_solidus, 'ro')
plot(T_fit, X_Cd_poly_fit, "k-")
plot(T_fit, X_Te_poly_fit, "k-")






%% Calc and Display Options  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

conditions.plot_log10_min = 12;   %%% set a threshold - if a quantity doesn't go above this at some temperature then dont bother plotting it
conditions.plot_log10_max = 20;   %%% set upper bound for plots too.  The real numbers are 10^x for these x's

% save to file options
% min cutoff concentration for significant concentrations
conditions.save_min = 1e12;  % set a minimum threshold for a defect or chargestate to be saved to the reduced output files (real number not log of it)

% define the material
conditions.stoich_flag = 'CdTe_native_only';
% conditions.stoich_flag = 'CdTe';

% save all the defects and chargestates or just the ones
% conditions.save_files_flag = 'All';
conditions.save_files_flag = 'Sig_Only';  % save only the ones above the save_min

% turn on Fermi-Dirac stats for bands or use Boltzmann stats
% conditions.Boltz_or_FD_flag = 'Boltz';
conditions.Boltz_or_FD_flag = 'FD';

% turn on or off site blocking statistics
conditions.site_blocking_flag= 'On_Infinite';  % site blocking for infinite crystal
% conditions.site_blocking_flag= 'On_Finite';  % site blocking for finite crystal
% conditions.site_blocking_flag= 'Off';

% turn on or off vibrational entropy.  Add 3kB per added atom, remove that for removed atoms.  G = H-TS so interstiitals should decrease G while vacacnies add it
% conditions.vib_ent_flag = '3kB';
% conditions.vib_ent_flag = 'Quantum';
conditions.vib_ent_flag = 'Off';

% Turn on or off T dependent band parameters
conditions.T_dep_bands_flag = 'On';
% conditions.T_dep_bands_flag = 'Off';

% give the relaxation energy for the STH
conditions.sth_flag = 0; % this turns off the STH's
% conditions.sth_flag = 1;  % this turns on the STH's
conditions.E_relax_sth1 = 0.0;  % relaxation energy for STH1's
conditions.E_relax_sth2 = 0.0;  % relaxation energy for STH2's

% set how to handle matrix chemical potentials
conditions.T_dep_matrix_mu_flag = 'On';      % using actual thermochemistry
% conditions.T_dep_matrix_mu_flag = 'Off';   % traditional way from most DFT papers - constant value vs T

% set how to handle defects with fixed concentrtations - either a onstant value for all T or different values vs T
% conditions.T_dep_fixed_defect_flag = 'On';
conditions.T_dep_fixed_defect_flag = 'Off';

% user can provide guesses for EF or Ef_mu_vec for each T_equilibrium that override automatic
% searching for grid searches.  For particleswarm, this guess is just added to the population of particles.
% conditions.Ef_mu_guesses_supplied_flag = 'On';
conditions.Ef_mu_guesses_supplied_flag = 'Off';
% conditions.Ef_mu_guesses = [];  % this makes the variable but the user needs to put the right values into it at some point below.
% conditions.Ef_mu_guesses = ones(11,1)*[1 -2];  % guesses for [Ef mufixed1 mufixed2...] one row per temperature

% choose the method of optimization used for finding solutions when elements are constrained - case 3 and 4
% conditions.search_method_flag = 'grid_fminsearch';
% conditions.search_method_flag = 'particleswarm_pattern';
conditions.search_method_flag = 'particleswarm_pattern_simplex';

conditions.search_method_quiet_flag = 'quiet';
% conditions.search_method_quiet_flag = 'verbose';

if strcmp(conditions.search_method_flag,'grid_fminsearch')
    %% set options for fminsearch for fixing elements concentrations
    % quicker set
    conditions.fixed_elements_fmin_MaxFunEvals = 1e5;
    conditions.fixed_elements_fmin_MaxIter = 1e5;
    conditions.fixed_elements_fmin_TolX = 1e-4;
    conditions.fixed_elements_fmin_TolFun = 1e-4;

    % finer set
    % conditions.fixed_elements_fmin_MaxFunEvals = 1e6;
    % conditions.fixed_elements_fmin_MaxIter = 1e6;
    % conditions.fixed_elements_fmin_TolX = 1e-5;
    % conditions.fixed_elements_fmin_TolFun = 1e-5;

elseif strcmp(conditions.search_method_flag,'particleswarm_pattern') || strcmp(conditions.search_method_flag,'particleswarm_pattern_simplex')
    % these control particle swarm search.  Assumption is that list of T's
    % goes from high to low.  For the first T, we do a more exhaustive
    % search for the solution.  Then, for temperatures 2-end we start out swarm2 at the prior solution and search within a band a few kBT wide.
    conditions.fixed_elements_swarm_iterlimit_base = 20;
    conditions.fixed_elements_swarm_stall_lim = 15;
    conditions.fixed_elements_swarm_display_int = 5;
    conditions.fixed_elements_swarm1_fine_factor = 3;
    conditions.fixed_elements_swarm1_min_size = 350;
    conditions.fixed_elements_swarm1_max_size = 15000;   %15000 seems to be ok upper limit on a laptop, 20000 is slow
    conditions.fixed_elements_swarm1_ObjectiveLimit = 1e10;  % quit if the best value is better than this value

    conditions.fixed_elements_swarm2_fine_factor = 3;
    conditions.fixed_elements_swarm2_search_band_kB = 2; % we take the solution from the prior temperature and search within a band +/- this*kB for the next one
    conditions.fixed_elements_swarm2_min_size = 300;
    conditions.fixed_elements_swarm2_max_size = 15000;
    conditions.fixed_elements_swarm2_ObjectiveLimit = 1e8;

    if strcmp(conditions.search_method_flag,'particleswarm_pattern_simplex')
        % these ones are for the final fminsearch after particleswarm is done.
        conditions.fixed_elements_fmin_MaxFunEvals = 1e4;
        conditions.fixed_elements_fmin_MaxIter = 7500;
        conditions.fixed_elements_fmin_TolX = 1e-5;
        conditions.fixed_elements_fmin_TolFun = 1e-5;
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Set any defects with fixed concentrations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conditions.fixed_defects = zeros(defects.num_defects,1);
conditions.fixed_defects_concentrations = zeros(defects.num_defects,1);

%% shallow dopants
conditions.Nd = 0;
conditions.Na = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% logic checks on fixed defetcs and sertting default to none fixed
% default is none, catch this.  Note that if elseif loops exit once one
% condition is satisfied

if ~isfield(conditions,'fixed_defects') || sum(isempty(conditions.fixed_defects))
    conditions.fixed_defects = zeros(defects.num_defects,1);
end

if strcmp(conditions.T_dep_fixed_defect_flag,'Off')
    conditions.fixed_defects_concentrations = conditions.fixed_defects_concentrations.*conditions.fixed_defects;          % if a defect is not fixed but by accident it has a target conc set, just set its target concentration to zero also, using the fixed defects vector as a mask
    if ~isfield(conditions,'fixed_defects_concentrations') || sum(isempty(conditions.fixed_defects_concentrations))>0
        conditions.fixed_defects_concentrations = zeros(defects.num_defects,1);
    end
elseif strcmp(conditions.T_dep_fixed_defect_flag,'On')
end

% applies to any fixed defect scenario
% find number of fixed defects and their indices
conditions.indices_of_fixed_defects = find(conditions.fixed_defects);
conditions.num_fixed_defects = numel(conditions.indices_of_fixed_defects);

if conditions.num_fixed_defects~=0
    for z=1:conditions.num_fixed_defects
        msg = strcat('WARNING: Fixed concentration used for...',defects.defect_names(conditions.indices_of_fixed_defects(z)),'...Be sure you meant to do this.');
        % disp(msg)
    end
elseif conditions.num_fixed_defects==0
    % disp('No Defects fixed.  Proceeding')
end

%%%  end section on fixed defects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Set up the elements to use in the model %%%%%%%%%%%%%
if ~isfield(conditions,'num_elements') || sum(isempty(conditions.num_elements))
    conditions.num_elements = size(defects.cs_dm,2);
end






% save copies of defects since we will be messing with it.
temp_defects = defects;






%% pressure
conditions.P_tot = 1;   % this will have no effect in this calculation as the pressures and mu's come into the mu later
conditions.P_units = 'atm';


%% Set T-dependnet properties like Eg and Vibent

% this has to come after the temperature vector is set
conditions.vibent_T0 = 200; % This is the characteristic To for the average phonon mode in the substance , as determined from the wo that makes the Debeye temperature come out to about 740 K based on the Cv(T) data found online.  To where x = hbar*w0/kBT = T0/T, wherein w0 is the mean phonon energy from the Debeye temperature determined from Cv(T).

% when T-independent values are used
conditions.TRef = 300;
conditions.EgRef = 1.5;   % this is for constant Eg(T)
conditions.EvRef = 0;
conditions.EcRef = conditions.EgRef;
conditions.NcRef = 1e18;
conditions.NvRef = 1e19;


% T-dependent Eg values
conditions.Eg0 = 1.5860;   % bandgap at 0K in eV.  Using Adrian's number from STEM EELS at PSU - 4.8 eV at 300 K.  So about 4.9 eV at 0K using
conditions.varshini_a = 5.9117e-4;  % can implement any model you want for Eg(T) as long as it gives a value for all T_equilibrium
conditions.varshini_b = 160;

% conditions.EcT_fraction = 0;
% conditions.EcT_fraction = 0.25;
% conditions.EcT_fraction = 0.5;
% conditions.EcT_fraction = 0.375;
% conditions.EcT_fraction = 0.50;     % what fraction of Eg(T)  happens in the CB?
% conditions.EcT_fraction = 0.75;
conditions.EcT_fraction = 0.80;   % From Intuon, it should be 80%
% conditions.EcT_fraction = 1;
conditions.EvT_fraction = 1-conditions.EcT_fraction;

% site densities in the lattice (num FU/unit cell / vol for unit cell)
N_Cd = 1.47e22;    %site 1
N_Te = 1.47e22;  % site 2
N_iCd = 1.47e22;  %site 3
N_iTe = 1.47e22;   %site 4
conditions.num_sites = [N_Cd; N_Te; N_iCd; N_iTe];   % note this needs to be a column vector, dont change it to a row vector

% % Calculate the numerical prefactor for each defect from temp_defects.degen_factor_config, temp_defects.degen_factor_elec, temp_defects.cs_num_each_site, and conditions.num_sites

%% temp_defects.cs_num_each_site tells you the number of primitive unit cells needed to form each defect (assuming the sites of the crystal, including distinct interstitials like the ia, ib, ... in b-ga2O3, are counted such that there is one per primitive unit cell
temp_defects.cs_site_prefactor = sort((ones(temp_defects.num_cs,1)*conditions.num_sites')./temp_defects.cs_num_each_site,2);      % calculates the max number of each chargestate that could be formed, given the availability of each site type in the crystal (or unit cells, or supercells, whatever the basis for counting).  Then sorts by the rows ascending, so the first column becomes the max number we could form.  Its done with a sort becasue you end up with Inf when you divide by sites you need zero of
temp_defects.cs_site_prefactor = temp_defects.cs_site_prefactor(:,1);  % toss out the other columns keeping only the first one, which should be the limiting one
temp_defects.cs_tot_prefactorfactor = temp_defects.cs_degen_factor_config .* temp_defects.cs_degen_factor_elec .* temp_defects.cs_site_prefactor;

clear N_Cd N_Te N_iCd N_iTe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



conditions.fixed_conc_flag = zeros(1,conditions.num_elements);
conditions.fixed_conc_values = zeros(1,conditions.num_elements);  % this will hold target values of elements
conditions.fixed_elements_mu_ranges = zeros(conditions.num_elements,2);   % create  matrix of right size for ALL the elements - only will pay attention to those which are fixed

conditions.mu_constant_flag = ones(1,conditions.num_elements);
conditions.mu_from_2nd_phases_flag = zeros(1,conditions.num_elements);   % set all to zero/false first then overwrite if needed
conditions.muT_equilibrium = -30*ones(1, temp_defects.num_elements); % set them all to extremely negative values then reset as needed
conditions.mu_set_flag = (conditions.mu_constant_flag + conditions.mu_from_2nd_phases_flag)>0;
neither_mu_nor_conc_flag = (conditions.mu_set_flag .* conditions.fixed_conc_flag);
both_mu_and_conc_flag = conditions.mu_set_flag + conditions.fixed_conc_flag ==0;

%% logic checks
mu_set_exclusively_by_one_method_check = ((conditions.mu_constant_flag + conditions.mu_from_2nd_phases_flag)-1) ~=0 ;
conditions.mu_set_flag = (conditions.mu_constant_flag + conditions.mu_from_2nd_phases_flag)>0;

if any(mu_set_exclusively_by_one_method_check)
    % disp(mu_set_exclusively_by_one_method_check)
    error('mu for each element has to be set either as a constant or by a phase boundary.  The line above has a 1 in the position of the element causing the problem')
end


% default is no elements fixed - catch if nothing entred.
if ~isfield(conditions,'fixed_conc_flag') || sum(isempty(conditions.fixed_conc_flag))
    conditions.fixed_conc_flag = zeros(temp_defects.num_elements,1);
elseif ~isfield(conditions,'fixed_conc_values') || sum(isempty(conditions.fixed_conc_values))
    conditions.fixed_conc_values = zeros(temp_defects.num_elements,1);
end
conditions.fixed_conc_values = conditions.fixed_conc_values.*conditions.fixed_conc_flag;
conditions.indices_of_fixed_elements = find(conditions.fixed_conc_flag);
conditions.num_fixed_elements = numel(conditions.indices_of_fixed_elements);
if conditions.num_fixed_elements~=0
    for z=1:conditions.num_fixed_elements
        msg = strcat('WARNING: Fixed concentration of element used for...',temp_defects.element_names(conditions.indices_of_fixed_elements(z)),'...Be sure you meant to do this.');
        % disp(msg)
    end
elseif conditions.num_fixed_elements==0
    % disp('No Eelements have fixed concentrations.  Proceeding')
end
clear msg impurity_mu_boundary_limits mu_set_exclusively_by_one_method_check max_mu_range start_kBT fine_factor
%%












%% calc X for Cd-rich first

% compute both sets of mu vs T, and the T values both can share
% note that each of these operations can change the range of T_fit becasue
% we only have liquidus data for some of the temperatures
[mu_T_Cd_rich, T_fit] = mu_CdTe_liquidus(T_fit, "Cd-rich");  % output is 2-column array [muCd muTe] and a 1 column vector of allowed T values
[mu_T_Te_rich, T_fit] = mu_CdTe_liquidus(T_fit, "Te-rich");  % output is 2-column array [muCd muTe] and a 1 column vector of allowed T values


N_temps = numel(T_fit);
output_holder = zeros(N_temps, 4 + defects.num_defects);   % T, total chi2, X_Cd_calc, X_Te_calc, total_Hform vector
lb = delta_E_lo * ones(1,defects.num_defects);    % this defines the bounds on the dE added to each defect's formation energy
ub = delta_E_hi * ones(1,defects.num_defects);

X_Cd_poly_fit = polyval(X_Cd_poly, T_fit);
X_Te_poly_fit = polyval(X_Te_poly, T_fit);


% swarm_opts = optimoptions('particleswarm','SwarmSize',swarmsize1,'ObjectiveLimit',Tloop_conditions.fixed_elements_swarm1_ObjectiveLimit,'MaxIterations',iter_lim,'MaxStallIterations',stall_lim,'HybridFcn', @patternsearch);
swarm_opts = optimoptions('particleswarm','Display','iter','HybridFcn', @patternsearch,'ObjectiveLimit',1e-15,'MaxStallIterations',10);
swarm_opts = optimoptions('particleswarm','HybridFcn', @patternsearch,'ObjectiveLimit',1e-15,'MaxStallIterations',10);




% loop through the list of T's
for i = 1:N_temps
    conditions.T_equilibrium = T_fit(i);

    conditions.kBT_equilibrium = conditions.kB * conditions.T_equilibrium;
    conditions.num_T_equilibrium = 1;
    conditions.NcT_equilibrium = conditions.NcRef*(conditions.T_equilibrium/conditions.TRef).^1.5;
    conditions.NvT_equilibrium = conditions.NvRef*(conditions.T_equilibrium/conditions.TRef).^1.5;
    conditions.delta_EgT = (conditions.varshini_a .* conditions.T_equilibrium.^2)./(conditions.T_equilibrium + conditions.varshini_b);  %delta is a positive number
    conditions.EgT_equilibrium = conditions.Eg0 - conditions.delta_EgT;
    conditions.EcT_equilibrium = conditions.Eg0*ones(size(conditions.EgT_equilibrium)) - conditions.EcT_fraction*conditions.delta_EgT;
    conditions.EvT_equilibrium = zeros(size(conditions.EgT_equilibrium)) + conditions.EvT_fraction*conditions.delta_EgT;


    fit_func = @(delta_H_form) fit_both(delta_H_form, X_Cd_poly_fit(i), X_Te_poly_fit(i), mu_T_Cd_rich(i,:), mu_T_Te_rich(i,:));

    [dH_optimal, GOF] = particleswarm(fit_func, defects.num_defects, lb, ub, swarm_opts);
    output_holder(i,1) = conditions.T_equilibrium;
    output_holder(i,2) = GOF;
    % set mu for the Cd_rich side, then compute X on Cd side
    conditions.muT_equilibrium = mu_T_Cd_rich(i,:);
    output_holder(i,3) = calc_X(dH_optimal);
    % set mu for the Te_rich side, then compute X on Te side
    conditions.muT_equilibrium = mu_T_Te_rich(i,:);
    output_holder(i,4) = calc_X(dH_optimal);
    output_holder(i,5:end) = dH_optimal;


end

figure(1)
plot(output_holder(:,1),output_holder(:,3),'go')
plot(output_holder(:,1),output_holder(:,4),'go')







    function [X_calc] = calc_X(delta_H_form)

        for m = 1:defects.num_defects
            temp_defects.cs_Eform(defects.cs_indices_lo(m):defects.cs_indices_hi(m)) = defects.cs_Eform(defects.cs_indices_lo(m):defects.cs_indices_hi(m)) + delta_H_form(m);
        end

        % run the calculaiton and calc stoichiometry
        [equilib_dark_sol] = defect_equilibrium_dark(conditions, temp_defects);
        [ratio_stoich, ~] = CdTe_stoich_native_only(equilib_dark_sol, conditions, temp_defects);
        X_calc = ratio_stoich(:,1);

    end


% function [chi2] = one_element_err(delta_H_form, X_target)
%     [X_calc] = calc_X(delta_H_form);
%     chi2 = (X_calc - X_target)^2;
% end



    function [chi2_tot] = fit_both(delta_H_form, X_Cd, X_Te, mu_T_Cd, mu_T_Te)

        % Cd part
        % set mu for the Cd_rich side, then compute err from Cd side
        conditions.muT_equilibrium = mu_T_Cd;
        for m = 1:defects.num_defects
            temp_defects.cs_Eform(defects.cs_indices_lo(m):defects.cs_indices_hi(m)) = defects.cs_Eform(defects.cs_indices_lo(m):defects.cs_indices_hi(m)) + delta_H_form(m);
        end
        [Cd_sol] = defect_equilibrium_dark(conditions, temp_defects);
        [ratio_stoich, ~] = CdTe_stoich_native_only(Cd_sol, conditions, temp_defects);
        chi2_Cd = ((ratio_stoich(:,1) - X_Cd)/X_Cd)^2;
        tot_defects_Cd = sum(Cd_sol.defects);

        % Te part
        % set mu for the Te_rich side, then compute err from Te side
        conditions.muT_equilibrium = mu_T_Te;
        for m = 1:defects.num_defects
            temp_defects.cs_Eform(defects.cs_indices_lo(m):defects.cs_indices_hi(m)) = defects.cs_Eform(defects.cs_indices_lo(m):defects.cs_indices_hi(m)) + delta_H_form(m);
        end
        [Te_sol] = defect_equilibrium_dark(conditions, temp_defects);
        [ratio_stoich, ~] = CdTe_stoich_native_only(Te_sol, conditions, temp_defects);
        chi2_Te = ((ratio_stoich(:,1) - X_Te)/X_Te)^2;
        tot_defects_Te = sum(Te_sol.defects);

        defects_penalty = (log10(tot_defects_Cd +  tot_defects_Te) - log10(5e19))^2;


        % sum the two errors plus penalties for having lots of defects
        % chi2_tot = (chi2_Cd + chi2_Te) + (tot_defects_Cd +  tot_defects_Te);
        chi2_tot = (chi2_Cd + chi2_Te) * defects_penalty;

    end 



end







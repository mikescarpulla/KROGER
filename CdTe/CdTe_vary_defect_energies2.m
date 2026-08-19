%% clear workspace and figures   %%%%%%%%%%%%
clear all
clc

%% set Physical Constants  %%%%%%%%%
conditions.q = 1.602176565e-19;
conditions.h = 6.62606957e-34;
conditions.kB = 8.6173324e-5;
conditions.mo = 9.1093837e-31;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Input the solidus line data - off stoichiometry boundaries for the CdTe solid phase
%% from Greenberg, J Crystal Growth 1996
%% Cd and Te compositions for the boundaries of the solidus region and the mu values for liquidus
%% T in Kelvin.  X from Greenberg originally in %, convert to 0<X<1.  Then to #/cm3 atoms
%% some extra points added to make both Cd and Te data span the whole T range and make them go to 50/50 at low T.  

T_Cd_rich = flip([600 697.2 819.6 864.1 871.3 900.8 939.7 969.4 992.0 1073 1123 1173 1223 1300]);    % list of measured points, changed to be decreasing T
pct_Cd_rich_solidus = flip([50 50.0002 50.0001 49.9990 49.9983 49.9990 49.9994 49.9991 49.9967 49.9966 49.9946 49.9957  49.9967 50]);
% pct_Cd_rich_errs = flip([ 0.0002 0.0001 0.0001 0.0006 0.0001 0.0002 0.0002 0.0003 0.0003 0.0003 0.0003 0.0003]);
T_Cd_rich = T_Cd_rich(1:end-1);  % delete the lowest-T point as we dont have liquidus data for it
pct_Cd_rich_solidus = pct_Cd_rich_solidus(1:end-1);
pct_Cd_rich_errs = pct_Cd_rich_errs(1:end-1);
X_Cd_rich_solidus = pct_Cd_rich_solidus'/100;   % this is now a fraction
X_Cd_rich_solidus_errs = pct_Cd_rich_errs'/100;    % this is now a fraction

%same for Te
T_Te_rich = flip([600 700 800 900 945.5 1016.6 1046.9 1073 1076.0 1094.1 1103.3 1121.0 1123 1173 1223 1243.1 1282.2 1285.5 1301.0 1312.0 1359.4 1361.8]);
pct_Te_rich_solidus = flip([50 50 50 50.0003000 50.0006 50.0036 50.0055 50.0071 50.0073 50.0092 50.0098 50.0138 50.0138 50.0135 50.0082 50.0087 50.0050 50.0047 50.0044 50.0038 50.0034 50.0013]);
% pct_Te_rich_errs =    flip([ 0.0001  0.0001  0.0004  0.0004  0.0001  0.0003  0.0001  0.0007  0.0007  0.0007  0.0007  0.0008  0.0010  0.0015  0.0030  0.0040  0.0006  0.0004]);
X_Te_rich_solidus = pct_Te_rich_solidus'/100;  % this is now in #/cm3
X_Te_rich_solidus_errs = pct_Te_rich_errs'/100;  % this is now in #/cm3

%% generate polynomial fits of the data that we can use for fitting calcs to.  
X_Cd_poly = polyfit(T_Cd_rich, X_Cd_rich_solidus, 4);
X_Te_poly = polyfit(T_Te_rich, X_Te_rich_solidus, 12);
X_Cd_poly_fit = polyval(X_Cd_poly, T_Cd_rich);
X_Te_poly_fit = polyval(X_Te_poly, T_Te_rich);

% plot the modified data and polynomial fits
figure(1)
clf
hold on
plot(T_Cd_rich, X_Cd_rich_solidus, 'bo')  %Greenberg data
plot(T_Te_rich, X_Te_rich_solidus, 'ro')
plot(T_Cd_rich, X_Cd_poly_fit, "k-")
plot(T_Te_rich, X_Te_poly_fit, "k-")


% choose lower and upper bounds
low_lim = 49.99;
hi_lim = 50.01;

%% vector of energies to add to the numbers in the defect database (can be negative or positive)
% these will be permuted across the different defects
% also how many error bars away should we allow candidates to be saved?
err_tol = 5;
% d_nrg = -2 : 0.25 : 2;
d_nrg = -3 : 0.5 : 3;
% d_nrg = flip(d_nrg);   % start from higher perturbations and work down
num_nrgs = numel(d_nrg);
conditions.save_fname = "CdTe_solidus_sweep_dHo.txt";



%%  Select the .mat file database of defect properties and load it %%%%%%%%

% [conditions.defect_db_fname,conditions.defect_db_pname] = uigetfile('*.mat','Select the defect database file you want to use');
conditions.defect_db_pname = "C:\Users\mscar\Documents\MATLAB\Defects";
conditions.defect_db_fname = "CdTe_Intuon_natives.mat";
% conditions.defect_db_fname = "CdTe_Intuon_vac_int_only.mat";








% figure out how many total scenarios (permuted formation energies) we want
% to check.  If there are too many, stop and tell user to try again.

num_scenarios = num_nrgs ^ defects.num_defects;

disp(num_scenarios)

input("confirm to start calc with this many scenarios? (hit enter)")




data_holder = zeros(num_scenarios, 5 + defects.num_cs);

% we will do this by using ndgrid to construct meshed arrays, then index
% into them using linear indexing.  what we will do is build a variabe command string.  Rather than making a
% fixed nested for loop, we can build the right number of permutations up front and then use the trick of linear indexing into multidimensional arrays to run the whole thing as one for loop.
% So build up a list of permuted scenarios up front, then just go through and run each one.



%% generate one array per defect that holds the perturbation of its formation energies for each permuted scenario
str1 = "";
str2 = "";
for j = 1:defects.num_defects
    str1 = str1 + "deltaE_" + int2str(j) + " ";
    str2= str2 + ", d_nrg";
end
str2 = char(str2); % strip off the first comma
str2 = string(str2(3:end));
text_command = "[" + str1 + "] = ndgrid(" + str2 + ");" ;
eval(text_command)

% at this point, we have generated a hyper-array (dimensions = num_defects) whcih holds permuted combinations of all the values in the d_nrg vector
% across all of the defects.  we can use linear indexing to pick out
% the d_nrg values for each defect and add those to the defects.dHo
% values for each of its chargestates for each scenario we run.
%
% The name of each array is deltaE_j where j runs from 1:num_defects








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






% save copies of defects and of conditions since we will be messing with
% them.  So below these lines everthing should be temp_
temp_defects = defects;
temp_conditions = conditions;

temp_conditions.P_tot = 1;   % this will have no effect in this calculation as the pressures and mu's come into the mu later
temp_conditions.P_units = 'atm';


% get the mu values vs T for both sides.  They dont change as we do all
% the scenarios
[mu_Cd_rich, T_out_Cd_rich] = mu_CdTe_liquidus(T_Cd_rich, "Cd-rich");  % output is 2-column array [muCd muTe]
[mu_Te_rich, T_out_Te_rich] = mu_CdTe_liquidus(T_Te_rich, "Te-rich");



%%%% Set T-dependnet properties like Eg and Vibent

% this has to come after the temperature vector is set
temp_conditions.vibent_T0 = 200; % This is the characteristic To for the average phonon mode in the substance , as determined from the wo that makes the Debeye temperature come out to about 740 K based on the Cv(T) data found online.  To where x = hbar*w0/kBT = T0/T, wherein w0 is the mean phonon energy from the Debeye temperature determined from Cv(T).

% when T-independent values are used
temp_conditions.TRef = 300;
temp_conditions.EgRef = 1.5;   % this is for constant Eg(T)
temp_conditions.EvRef = 0;
temp_conditions.EcRef = temp_conditions.EgRef;
temp_conditions.NcRef = 1e18;
temp_conditions.NvRef = 1e19;


% T-dependent Eg values
temp_conditions.Eg0 = 1.5860;   % bandgap at 0K in eV.  Using Adrian's number from STEM EELS at PSU - 4.8 eV at 300 K.  So about 4.9 eV at 0K using
temp_conditions.varshini_a = 5.9117e-4;  % can implement any model you want for Eg(T) as long as it gives a value for all T_equilibrium
temp_conditions.varshini_b = 160;

% temp_conditions.EcT_fraction = 0;
% temp_conditions.EcT_fraction = 0.25;
% temp_conditions.EcT_fraction = 0.5;
% temp_conditions.EcT_fraction = 0.375;
% temp_conditions.EcT_fraction = 0.50;     % what fraction of Eg(T)  happens in the CB?
% temp_conditions.EcT_fraction = 0.75;
temp_conditions.EcT_fraction = 0.80;   % From Intuon, it should be 80%
% temp_conditions.EcT_fraction = 1;
temp_conditions.EvT_fraction = 1-temp_conditions.EcT_fraction;

% site densities in the lattice (num FU/unit cell / vol for unit cell)
N_Cd = 1.47e22;    %site 1
N_Te = 1.47e22;  % site 2
N_iCd = 1.47e22;  %site 3
N_iTe = 1.47e22;   %site 4
temp_conditions.num_sites = [N_Cd; N_Te; N_iCd; N_iTe];   % note this needs to be a column vector, dont change it to a row vector

% % Calculate the numerical prefactor for each defect from temp_defects.degen_factor_config, temp_defects.degen_factor_elec, temp_defects.cs_num_each_site, and temp_conditions.num_sites

%% temp_defects.cs_num_each_site tells you the number of primitive unit cells needed to form each defect (assuming the sites of the crystal, including distinct interstitials like the ia, ib, ... in b-ga2O3, are counted such that there is one per primitive unit cell
temp_defects.cs_site_prefactor = sort((ones(temp_defects.num_cs,1)*temp_conditions.num_sites')./temp_defects.cs_num_each_site,2);      % calculates the max number of each chargestate that could be formed, given the availability of each site type in the crystal (or unit cells, or supercells, whatever the basis for counting).  Then sorts by the rows ascending, so the first column becomes the max number we could form.  Its done with a sort becasue you end up with Inf when you divide by sites you need zero of
temp_defects.cs_site_prefactor = temp_defects.cs_site_prefactor(:,1);  % toss out the other columns keeping only the first one, which should be the limiting one
temp_defects.cs_tot_prefactorfactor = temp_defects.cs_degen_factor_config .* temp_defects.cs_degen_factor_elec .* temp_defects.cs_site_prefactor;

clear N_Cd N_Te N_iCd N_iTe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% compute 2 sets of T dependent Egap numbers
Nc_Cd_rich = temp_conditions.NcRef*(T_out_Cd_rich/temp_conditions.TRef).^1.5;
Nv_Cd_rich = temp_conditions.NvRef*(T_out_Cd_rich/temp_conditions.TRef).^1.5;
delta_Eg_Cd_rich = (temp_conditions.varshini_a .* T_out_Cd_rich.^2)./(T_out_Cd_rich + temp_conditions.varshini_b);  %delta is a positive number
Eg_Cd_rich = temp_conditions.Eg0 - delta_Eg_Cd_rich;
Ec_Cd_rich = temp_conditions.Eg0*ones(size(Eg_Cd_rich)) - temp_conditions.EcT_fraction*delta_Eg_Cd_rich;
Ev_Cd_rich = zeros(size(Eg_Cd_rich)) + temp_conditions.EvT_fraction*delta_Eg_Cd_rich;

Nc_Te_rich = temp_conditions.NcRef*(T_out_Te_rich/temp_conditions.TRef).^1.5;
Nv_Te_rich = temp_conditions.NvRef*(T_out_Te_rich/temp_conditions.TRef).^1.5;
delta_Eg_Te_rich = (temp_conditions.varshini_a .* T_out_Te_rich.^2)./(T_out_Te_rich + temp_conditions.varshini_b);  %delta is a positive number
Eg_Te_rich = temp_conditions.Eg0 - delta_Eg_Te_rich;
Ec_Te_rich = temp_conditions.Eg0*ones(size(Eg_Te_rich)) - temp_conditions.EcT_fraction*delta_Eg_Te_rich;
Ev_Te_rich = zeros(size(Eg_Te_rich)) + temp_conditions.EvT_fraction*delta_Eg_Te_rich;






% now start the calculation.  We are still inside the i=1,2 loop for Cd- or
% Te- rich cases

% figure(1)
% clf
% hold on
% chi2 = plot(0,0,'bo');

counter = 0;

for l = 1:num_scenarios  % we can index in to a 3x3 array 1 = (1,1,1), 2=(1,1,2),.. 4 = (1,2,1), 5=(1,2,2) etc - we jsut do it in N dimensions now.  It doesnt matter the order we check them as longas we get the corresponding dHo from all of the defects systematically
    % disp(int2str(l))

    for m = 1:defects.num_defects
        for n = defects.cs_indices_lo(m) : defects.cs_indices_hi(m)   % index over the chargestates of the jth defect
            temp_defects.cs_Eform(n) = defects.cs_Eform(n) + eval("deltaE_" + int2str(m) + "(" + int2str(l) + ")");  %find the dHo of the kth chargestate and add the jth delta_energ to it. Store the result in a temprary defects varaible that changes for each scenario
        end
    end
    % at this point we havechanged all the Eform values for this scenario

    temp_conditions.fixed_conc_flag = zeros(1,temp_conditions.num_elements);
    temp_conditions.fixed_conc_values = zeros(1,temp_conditions.num_elements);  % this will hold target values of elements
    temp_conditions.fixed_elements_mu_ranges = zeros(temp_conditions.num_elements,2);   % create  matrix of right size for ALL the elements - only will pay attention to those which are fixed



    %% do Cd and Te rich cases

    for i=1:2  % loop over Cd rich i=1 or Te rich i=2 cases

        if i==1

            temp_conditions.T_equilibrium = T_out_Cd_rich;
            temp_conditions.num_T_equilibrium = numel(temp_conditions.T_equilibrium);
            temp_conditions.kBT_equilibrium = temp_conditions.kB * temp_conditions.T_equilibrium;

            temp_conditions.mu_constant_flag = ones(1,temp_conditions.num_elements);
            temp_conditions.mu_from_2nd_phases_flag = zeros(1,temp_conditions.num_elements);   % set all to zero/false first then overwrite if needed
            temp_conditions.muT_equilibrium = -30*ones(temp_conditions.num_T_equilibrium, temp_defects.num_elements); % set them all to extremely negative values then reset as needed
            temp_conditions.mu_set_flag = (temp_conditions.mu_constant_flag + temp_conditions.mu_from_2nd_phases_flag)>0;
            neither_mu_nor_conc_flag = (temp_conditions.mu_set_flag .* temp_conditions.fixed_conc_flag);
            both_mu_and_conc_flag = temp_conditions.mu_set_flag + temp_conditions.fixed_conc_flag ==0;

            temp_conditions.muT_equilibrium = mu_Cd_rich;

            temp_conditions.NcT_equilibrium = Nc_Cd_rich;
            temp_conditions.NvT_equilibrium = Nv_Cd_rich;
            delta_EgT_equilibrium = delta_Eg_Cd_rich;
            temp_conditions.EgT_equilibrium = Eg_Cd_rich;
            temp_conditions.EcT_equilibrium = Ec_Cd_rich;
            temp_conditions.EvT_equilibrium = Ev_Cd_rich;


        elseif i==2

            temp_conditions.T_equilibrium = T_out_Te_rich;
            temp_conditions.num_T_equilibrium = numel(temp_conditions.T_equilibrium);
            temp_conditions.kBT_equilibrium = temp_conditions.kB * temp_conditions.T_equilibrium;

            temp_conditions.mu_constant_flag = ones(1,temp_conditions.num_elements);
            temp_conditions.mu_from_2nd_phases_flag = zeros(1,temp_conditions.num_elements);   % set all to zero/false first then overwrite if needed
            temp_conditions.muT_equilibrium = -30*ones(temp_conditions.num_T_equilibrium, temp_defects.num_elements); % set them all to extremely negative values then reset as needed
            temp_conditions.mu_set_flag = (temp_conditions.mu_constant_flag + temp_conditions.mu_from_2nd_phases_flag)>0;
            neither_mu_nor_conc_flag = (temp_conditions.mu_set_flag .* temp_conditions.fixed_conc_flag);
            both_mu_and_conc_flag = temp_conditions.mu_set_flag + temp_conditions.fixed_conc_flag ==0;

            temp_conditions.muT_equilibrium = mu_Te_rich;

            temp_conditions.NcT_equilibrium = Nc_Te_rich;
            temp_conditions.NvT_equilibrium = Nv_Te_rich;
            delta_EgT_equilibrium = delta_Eg_Te_rich;
            temp_conditions.EgT_equilibrium = Eg_Te_rich;
            temp_conditions.EcT_equilibrium = Ec_Te_rich;
            temp_conditions.EvT_equilibrium = Ev_Te_rich;

        end
        %% now we have the T vectors and mu arrays for the Cd or Te rich sides





        %% logic checks
        mu_set_exclusively_by_one_method_check = ((temp_conditions.mu_constant_flag + temp_conditions.mu_from_2nd_phases_flag)-1) ~=0 ;
        temp_conditions.mu_set_flag = (temp_conditions.mu_constant_flag + temp_conditions.mu_from_2nd_phases_flag)>0;

        if any(mu_set_exclusively_by_one_method_check)
            % disp(mu_set_exclusively_by_one_method_check)
            error('mu for each element has to be set either as a constant or by a phase boundary.  The line above has a 1 in the position of the element causing the problem')
        end


        % default is no elements fixed - catch if nothing entred.
        if ~isfield(conditions,'fixed_conc_flag') || sum(isempty(temp_conditions.fixed_conc_flag))
            temp_conditions.fixed_conc_flag = zeros(temp_defects.num_elements,1);
        elseif ~isfield(conditions,'fixed_conc_values') || sum(isempty(temp_conditions.fixed_conc_values))
            temp_conditions.fixed_conc_values = zeros(temp_defects.num_elements,1);
        end
        temp_conditions.fixed_conc_values = temp_conditions.fixed_conc_values.*temp_conditions.fixed_conc_flag;
        temp_conditions.indices_of_fixed_elements = find(temp_conditions.fixed_conc_flag);
        temp_conditions.num_fixed_elements = numel(temp_conditions.indices_of_fixed_elements);
        if temp_conditions.num_fixed_elements~=0
            for z=1:temp_conditions.num_fixed_elements
                msg = strcat('WARNING: Fixed concentration of element used for...',temp_defects.element_names(temp_conditions.indices_of_fixed_elements(z)),'...Be sure you meant to do this.');
                % disp(msg)
            end
        elseif temp_conditions.num_fixed_elements==0
            % disp('No Eelements have fixed concentrations.  Proceeding')
        end
        clear msg impurity_mu_boundary_limits mu_set_exclusively_by_one_method_check max_mu_range start_kBT fine_factor



        % run the calculaiton at all the temperatures for Cd or Te rich
        % (i=1 or i=2 cases)
        [equilib_dark_sol] = defect_equilibrium_dark(temp_conditions, temp_defects);

        % calc stoichiometry
        [ratio_stoich, ~] = CdTe_stoich_native_only(equilib_dark_sol, temp_conditions, temp_defects);
        X_calc = ratio_stoich(:,1);



        % check if this entire scenario is within N times the error from Goldberg for all of the temperatures

        majority = ceil(numel(X_calc)/2);

        if i==1 %Cd-rich
            max_X_lim = X_Cd_rich_solidus + err_tol * X_Cd_rich_solidus_errs;
            min_X_lim = X_Cd_rich_solidus - err_tol * X_Cd_rich_solidus_errs;
            % Cd_rich_delta_X = X_calc-0.5;
            Cd_rich_close_flag = all([sum(X_calc <= max_X_lim)>=majority sum(X_calc >= min_X_lim)>=majority]);
            Cd_rich_chi2 = 1/(temp_conditions.num_T_equilibrium-1) * sum(  ((X_calc - X_Cd_rich_solidus)./X_Cd_rich_solidus_errs).^2 );
        elseif i==2 %Te-rich
            max_X_lim = X_Te_rich_solidus + err_tol * X_Te_rich_solidus_errs;
            min_X_lim = X_Te_rich_solidus - err_tol * X_Te_rich_solidus_errs;
            Te_rich_close_flag = all([sum(X_calc <= max_X_lim)>=majority sum(X_calc >= min_X_lim)>=majority]);
            Te_rich_chi2 = 1/(temp_conditions.num_T_equilibrium-1) * sum(  ((X_calc - X_Te_rich_solidus)./X_Te_rich_solidus_errs).^2 );
        end



    end   %end loop over i

    % % update the plot
    % chi2.YData = Cd_rich_chi2 + Te_rich_chi2;
    % chi2.XData = l;
    % drawnow limitrate;





    % only save if the Cd and Te sides are within the limits
    if Cd_rich_close_flag && Te_rich_close_flag
        % disp(["scenario ", num2str(l)])
        counter = counter + 1;
        total_chi2 = Cd_rich_chi2 + Te_rich_chi2;
        % disp(["successful scenario ",int2str(counter)])
        data_holder(counter,:) = [counter l Cd_rich_chi2 Te_rich_chi2 total_chi2 temp_defects.cs_Eform'];
    end

end  %end l loop


% trim off all of the rows of all zeros at the end of the pre-allotted data holder
data_holder = data_holder(1:counter,:);
data_holder = sortrows(data_holder,5);



% write output file
headers = ["successful scenario" "absolute scenario" "Cd_rich_chi2" "Te_rich_chi2" "total_chi2" temp_defects.cs_names'];
output_cell = [cellstr(headers); num2cell(data_holder)];
writecell(output_cell,strcat(temp_conditions.save_pname, temp_conditions.save_fname),'Delimiter','tab','WriteMode','overwrite');
disp('wrote summaries of each scenario to text file')


% Plot the best scenarios
figure(1)
clf
hold on
plot(X_Cd_rich_solidus, T_Cd_rich, 'bo')  %Greenberg data
plot(X_Te_rich_solidus, T_Te_rich, 'ro')



% sort the results by the total chi2.  default direction is ascending
% values
topX = 20;
top_results = data_holder(1:topX, :);  % take only the best topX fits

Tapprox = 1360:-60:820;  % approx range, this will be updated by the net two lines.
[mu_Cd_rich, T_out_Cd_rich] = mu_CdTe_liquidus(Tapprox, "Cd-rich");  % output is 2-column array [muCd muTe]
[mu_Te_rich, T_out_Te_rich] = mu_CdTe_liquidus(Tapprox, "Te-rich");

 


for k=1:topX
    scenario = top_results(k,1)
    temp_defects.cs_Eform = (top_results(k, 6:end))';

    for j=1:2
        if j==1
            temp_conditions.T_equilibrium = T_out_Cd_rich;
            temp_conditions.num_T_equilibrium = numel(temp_conditions.T_equilibrium);
            temp_conditions.kBT_equilibrium = temp_conditions.kB * temp_conditions.T_equilibrium;

            temp_conditions.muT_equilibrium = mu_Cd_rich;
            temp_conditions.NcT_equilibrium = temp_conditions.NcRef*(T_out_Cd_rich/temp_conditions.TRef).^1.5;
            temp_conditions.NvT_equilibrium = temp_conditions.NvRef*(T_out_Cd_rich/temp_conditions.TRef).^1.5;
            temp_conditions.delta_EgT_equilibrium = (temp_conditions.varshini_a .* T_out_Cd_rich.^2)./(T_out_Cd_rich + temp_conditions.varshini_b);  %delta is a positive number
            temp_conditions.EgT_equilibrium = temp_conditions.Eg0 - temp_conditions.delta_EgT_equilibrium;
            temp_conditions.EcT_equilibrium = temp_conditions.Eg0*ones(size(temp_conditions.delta_EgT_equilibrium)) - temp_conditions.EcT_fraction*temp_conditions.delta_EgT_equilibrium;
            temp_conditions.EvT_equilibrium = zeros(size(temp_conditions.delta_EgT_equilibrium)) + temp_conditions.EvT_fraction*temp_conditions.delta_EgT_equilibrium;

            temp_conditions.mu_constant_flag = ones(1,temp_conditions.num_elements);
            temp_conditions.mu_from_2nd_phases_flag = zeros(1,temp_conditions.num_elements);   % set all to zero/false first then overwrite if needed
            temp_conditions.muT_equilibrium = -30*ones(temp_conditions.num_T_equilibrium, temp_defects.num_elements); % set them all to extremely negative values then reset as needed
            temp_conditions.mu_set_flag = (temp_conditions.mu_constant_flag + temp_conditions.mu_from_2nd_phases_flag)>0;
            neither_mu_nor_conc_flag = (temp_conditions.mu_set_flag .* temp_conditions.fixed_conc_flag);
            both_mu_and_conc_flag = temp_conditions.mu_set_flag + temp_conditions.fixed_conc_flag ==0;

            [Cd_rich_sol] = defect_equilibrium_dark(temp_conditions, temp_defects);
            [ratio_stoich, ~] = CdTe_stoich_native_only(Cd_rich_sol, temp_conditions, temp_defects);
            X_Cd_rich = ratio_stoich(:,1);
            plot(X_Cd_rich, T_out_Cd_rich)

        elseif j==2
            temp_conditions.T_equilibrium = T_out_Te_rich;
            temp_conditions.num_T_equilibrium = numel(temp_conditions.T_equilibrium);
            temp_conditions.kBT_equilibrium = temp_conditions.kB * temp_conditions.T_equilibrium;

            temp_conditions.muT_equilibrium = mu_Te_rich;
            temp_conditions.NcT_equilibrium = temp_conditions.NcRef*(T_out_Te_rich/temp_conditions.TRef).^1.5;
            temp_conditions.NvT_equilibrium = temp_conditions.NvRef*(T_out_Te_rich/temp_conditions.TRef).^1.5;
            temp_conditions.delta_EgT_equilibrium = (temp_conditions.varshini_a .* T_out_Te_rich.^2)./(T_out_Te_rich + temp_conditions.varshini_b);  %delta is a positive number
            temp_conditions.EgT_equilibrium = temp_conditions.Eg0 - temp_conditions.delta_EgT_equilibrium;
            temp_conditions.EcT_equilibrium = temp_conditions.Eg0*ones(size(temp_conditions.delta_EgT_equilibrium)) - temp_conditions.EcT_fraction*temp_conditions.delta_EgT_equilibrium;
            temp_conditions.EvT_equilibrium = zeros(size(temp_conditions.delta_EgT_equilibrium)) + temp_conditions.EvT_fraction*temp_conditions.delta_EgT_equilibrium;

            temp_conditions.mu_constant_flag = ones(1,temp_conditions.num_elements);
            temp_conditions.mu_from_2nd_phases_flag = zeros(1,temp_conditions.num_elements);   % set all to zero/false first then overwrite if needed
            temp_conditions.muT_equilibrium = -30*ones(temp_conditions.num_T_equilibrium, temp_defects.num_elements); % set them all to extremely negative values then reset as needed
            temp_conditions.mu_set_flag = (temp_conditions.mu_constant_flag + temp_conditions.mu_from_2nd_phases_flag)>0;
            neither_mu_nor_conc_flag = (temp_conditions.mu_set_flag .* temp_conditions.fixed_conc_flag);
            both_mu_and_conc_flag = temp_conditions.mu_set_flag + temp_conditions.fixed_conc_flag ==0;

            [Te_rich_sol] = defect_equilibrium_dark(temp_conditions, temp_defects);
            [ratio_stoich, ~] = CdTe_stoich_native_only(Te_rich_sol, temp_conditions, temp_defects);
            X_Te_rich = ratio_stoich(:,1);
            plot(X_Te_rich, T_out_Te_rich)
        end
    end
end



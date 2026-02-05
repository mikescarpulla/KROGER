% T in Kelvin, kBT in eV, and the DOS factors which are set so that
% the integrals for n,p can be run in eV units and come out to #/cm3
% the example here is Ga2O3 with defects computed by Joel Varley
% convetion for sites: Ga1 Ga2 O1 O2 O3 interstitial ...
% convention for chemical potentials [Ga O Si H Fe Sn ...]
% Jan 2023 - added ability to freeze defect concentrations
%feb 2023 - added site blocking
% 2/14/2023 - added provisions for variable pO2 - can make brower diagram plots next
% 3/23 - 1/2024 added ability to freeze elements, major code overhauls -
% Megan Stephans and MS.
% 3/28/24 - renamed variables in defects to cs_ to denote they have to do
% with the chargestates.  Implemented ability for calculation of prefactors for each
% chargestate using electronic degeneracy, config degeneracy, and number of
% primitive cells needed (complexes including multiple Sn_Ga defects would
% eat up multiple primitive cells to form each one).
% 4/3/24 - replaced grid search + fminsearch with particle swarm +
% patternsearch, tuned default parameters.  Fixed some issues with saving
% and the summed error function.  Fixed saving files.  Fixed bugs.  Tuned
% parameters for swarm optimization in high dimensions without running out
% of memory
% 4/8/24 - added reduced files to save for defects above a user defined threshold conc for all T's .  Cleaned up plotting.  Added ability to plot fractions of fixed elements in each defect.
% 8/12/24 - added self trapped holes as another flavor of holes.  Added
% quantum dSvib based on number of atoms changed out

%% To Do items: make the site blocking work by site differently not just total

%% clear workspace and figures   %%%%%%%%%%%%
clear all
clc

conditions.figures_flag = 'On';
% conditions.figures_flag = 'Off';

% if strcmp(conditions.figures_flag,'On')
%     for i=1:7
%         figure(i)
%         clf
%         hold on
%     end
%     clear i
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% set Physical Constants  %%%%%%%%%
conditions.q = 1.602176565e-19;
conditions.h = 6.62606957e-34;
conditions.kB = 8.6173324e-5;
conditions.mo = 9.1093837e-31;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set problem definition and calc details by editing this section %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% future improvement - make this into a modular ingestion where you can
%%% mix and match different subsets of all the defects and it will number
%%% the chargestates automatically


%%  Select the .mat file database of defect properties and load it %%%%%%%%
%%% note: defects and chargestates listed down, one row per.  Elements
%%% listed as columns in variables like defects.cs_dm.  These conventions must
%%% be obeyed so that the matrix and vector multiplicaitons used in calcs
%%% work out right - e.g. when mu vectors and dm vectors are multiplied to
%%% get formation energies for each defect.

conditions.defect_db_fname = [];
conditions.defect_db_pname = [];
% [conditions.defect_db_fname,conditions.defect_db_pname] = uigetfile('*.mat','Select the defect database file you want to use');

if ~ischar(conditions.defect_db_fname) || isempty(conditions.defect_db_fname)   %if the user doesnt pick one in the UIgetfile, do it manually here
    % conditions.defect_db_fname = 'Ga2O3_Varley_all_defects_new_072924.mat';
    conditions.defect_db_fname = 'Ga2O3_Varley_all_defects_new_092724.mat';
end

if ~ischar(conditions.defect_db_pname) || isempty(conditions.defect_db_pname)
    conditions.defect_db_pname = strcat(pwd,'\');  % the pwd command leaves off the last \ so have to add it.
end

load(strcat(conditions.defect_db_pname,conditions.defect_db_fname))
defect_db_name_root = conditions.defect_db_fname(1:(size(conditions.defect_db_fname,2) - 4)    );  % this takes the filename and strips the . and extension - keep it until savenames established

if numel(fieldnames(defects)) < 18
    error('Defects database file seems to be missing some required parts - carefully check it and try again')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% 
% %%%%%%%%%%%%%%%%%% Manually increase dHo for Vga and complexes - numbers based on Varley 031424 database
% for ii = [38:243 252:559 568:605 632:675]
%     % defects.cs_dHo(ii) = defects.cs_dHo(ii) + 0.5;
%     defects.cs_dHo(ii) = defects.cs_dHo(ii) - 0.5;
%     % defects.cs_dHo(ii) = defects.cs_dHo(ii) - 1;
% end
% clear ii
% %%%%%%%%%%%%%%%%





%%%%% choose folder and base filename to save calculation run outputs
[save_fname_root,conditions.save_pname] = uiputfile('*.*','Select or type the filename and location for saving results');
save_fname_root_length = size(save_fname_root,2);

if ~ischar(save_fname_root)  %if user exits the UI window with no input
    save_fname_root = strcat(defect_db_name_root,strcat('_',strrep(char(datetime),' ','_')));  % default is to use the defect database name but append date and time to it
    conditions.save_pname = conditions.defect_db_pname;  %default is to use the same directory with the defect database in it
elseif strcmp(save_fname_root(save_fname_root_length-3),'.')
    save_fname_root = save_fname_root(1:(save_fname_root_length-4) );   % this takes the filename and strips the extension.  If you type a name without extension, the uiputfile adds .dat to it
    conditions.save_pname = conditions.defect_db_pname;  %default is to use the same directory with the defect database in it
end

conditions.conditions_save_fname = strcat(save_fname_root,'_calculation_conditions');
conditions.equilib_save_fname = strcat(save_fname_root,'_equilib_dark_sol');
conditions.fullquench_save_fname = strcat(save_fname_root,'_fullquench_dark_sol');
conditions.elements_save_fname = strcat(save_fname_root,'_element_totals');
conditions.grouping_save_fname = strcat(save_fname_root,'_grouped_defects');
clear save_fname_root save_fname_root_length defect_set_flag defect_db_name_root
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Calc and Display Options  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

conditions.plot_log10_min = 14;   %%% set a threshold - if a quantity doesn't go above this at some temperature then dont bother plotting it
conditions.plot_log10_max = 19;   %%% set upper bound for plots too.  The real numbers are 10^x for these x's

% save to file options
% min cutoff concentration for significant concentrations
conditions.save_min = 1e14;  % set a minimum threshold for a defect or chargestate to be saved to the reduced output files (real number not log of it)

% conditions.stoich_flag = 'CdTe';
conditions.stoich_flag = 'Ga2O3';

% save all the defects and chargestates or just the ones
% conditions.save_files_flag = 'All';
conditions.save_files_flag = 'Sig_Only';  % save only the ones above the save_min

% conditions.defect_group_flag = 'On';
conditions.defect_group_flag = 'Off';
if strcmp(conditions.defect_group_flag,'On')
    % conditions.defect_grouping_fname = 'CdTe_grouping_defects.mat';
    conditions.defect_grouping_fname = 'Ga2O3_grouping_defects.mat';
    conditions.defect_grouping_pname = strcat(pwd,'\');
end

% turn on Fermi-Dirac stats for bands or use Boltzmann stats
% conditions.Boltz_or_FD_flag = 'Boltz';
conditions.Boltz_or_FD_flag = 'FD';

% turn on or off site blocking statistics
conditions.site_blocking_flag= 'On_Infinite';  % site blocking for infinite crystal
% conditions.site_blocking_flag= 'On_Finite';  % site blocking for finite crystal
% conditions.site_blocking_flag= 'Off';

% turn on or off vibrational entropy.  Add 3kB per added atom, remove that for removed atoms.  G = H-TS so interstiitals should decrease G while vacacnies add it
% conditions.vib_ent_flag = '3kB';
conditions.vib_ent_flag = 'Quantum';
% conditions.vib_ent_flag = 'Off';

% Turn on or off T dependent band parameters
conditions.T_dep_bands_flag = 'On';
% conditions.T_dep_bands_flag = 'Off';

% give the relaxation energy for the STH
% conditions.sth_flag = 0; % this turns off the STH's
conditions.sth_flag = 1;  % this turns on the STH's
conditions.E_relax_sth1 = 0.52;  % relaxation energy for STH1's
conditions.E_relax_sth2 = 0.53;  % relaxation energy for STH2's


% set how to handle matrix chemical potentials
conditions.T_dep_matrix_mu_flag = 'On';      % using actual thermochemistry
% conditions.T_dep_matrix_mu_flag = 'Off';   % traditional way from most DFT papers - constant value vs T

% set how to handle defects with fixed concentrtations - either a
% constant value for all T or different values vs T
% conditions.T_dep_fixed_defect_flag = 'On';
conditions.T_dep_fixed_defect_flag = 'Off';

% paraequilibrium is when we assume that some of the defects like
% interstitials will equilibrate like electrons in quenching.
% conditions.interstitial_equilibrate_at_Tquench_flag = 'On';
% conditions.interstitial_equilibrate_at_Tquench_flag = 'Off';

% user can provide guesses for EF or Ef_mu_vec for each T_equilibrium that override automatic
% searching for grid searches.  For particleswarm, this guess is just added to the population of particles.  
% conditions.Ef_mu_guesses_supplied_flag = 'On';
conditions.Ef_mu_guesses_supplied_flag = 'Off';
% conditions.Ef_mu_guesses = [];  % this makes the variable but the user needs to put the right values into it at some point below.   
% conditions.Ef_mu_guesses = ones(11,1)*[1 -2];  % guesses for [Ef mufixed1 mufixed2...] one row per temperature

% adjust how fine the steps are for brute force search over Ef for quenching.  spacing = kBT/this.  kBT/2 is too coarse a lot of the time
conditions.fullquench_EF_search_step_divisor = 5;

% choose the method of optimization used for finding solutions when elements are constrained - case 3 and 4
% conditions.search_method_flag = 'grid_fminsearch';
% conditions.search_method_flag = 'particleswarm_pattern';
conditions.search_method_flag = 'particleswarm_pattern_simplex';

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



%% set up batch calculation permuting 4 dopants and 2 values of f
%
% for index = 1:8
%
%     if index==1 || index==5
%         save_fname_root = 'Sn_Kuramata_2016';
%     elseif index==2 || index==6
%         save_fname_root = 'Mg_doped';
%     elseif index==3 || index==7
%         save_fname_root = 'Fe_doped';
%     elseif index==4 || index==8
%         save_fname_root = 'UID_Kuramata_2016';
%     end
%
%     if index<=4
%         conditions.EcT_fraction = 0.75;     % what fraction of Eg change vs T happens in the CB?
%         save_fname_root = strcat(save_fname_root,'_f75');
%     elseif index>=5
%         conditions.EcT_fraction = 0.50;     % what fraction of Eg change vs T happens in the CB?
%         save_fname_root = strcat(save_fname_root,'_f50');
%     end
%
%     conditions_save_fname = strcat(save_fname_root,'_calculation_conditions');
%     equilib_save_fname = strcat(save_fname_root,'_equilib_dark_sol');
%     fullquench_save_fname = strcat(save_fname_root,'_fullquench_dark_sol');
%
%
%     if strcmp(conditions.figures_flag,'On')
%         for i=1:7
%             figure(i)
%             clf
%             hold on
%         end
%         clear i
%         conditions.plot_log10_min = 14;   %%% set a threshold - if a quantity doesn't go above this at some temperature then dont bother plotting it
%         conditions.plot_log10_max = 20;   %%% set upper bound for plots too
%     end
%
%



%% set Temperatures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vector of T values at which to calculate for equilibrium
% start at the highest T you want and then work down in order to fully
% exploit the adaptive guessing.  Also, for fixed elements the major time
% is spent finding the first guess, the optimization is fast by comparison once you fix 2 or more elements.
% So it's worth doing a fine T grid for long calcuations so you dont have
% to ever do it over.

% conditions.T_equilibrium = 1500:-100:700;  % row vector of T in K at which to compute the equilirium
% conditions.T_equilibrium = 2100:-100:500;
% conditions.T_equilibrium = 2100:-100:1500;
% conditions.T_equilibrium = 1200:-50:800;
% conditions.T_equilibrium = 2100:-100:2000;
% conditions.T_equilibrium = 10:10:300;
% conditions.T_equilibrium = 2100:-50:1400;
% conditions.T_equilibrium = 500:50:2100;
% conditions.T_equilibrium = 500:100:2100;
% conditions.T_equilibrium = 300:100:2100;


conditions.T_equilibrium = 1000 + 273;   % use this for tasks like doping or mu sweeps done at one temperature
% conditions.T_equilibrium = 600 + 273;
% conditions.T_equilibrium = 1600 + 273;
% conditions.T_equilibrium = 1200 + 273;
% conditions.T_equilibrium = 1100 + 273;
% conditions.T_equilibrium = 1600 + 273;

conditions.num_T_equilibrium = numel(conditions.T_equilibrium);
conditions.kBT_equilibrium = conditions.kB * conditions.T_equilibrium;

% set final T for full quenching from the T_equilibrium values above
conditions.T_fullquench = 300;   % scalar temperature at which to compute the defect concentrations after quenching from T_equilibrium
conditions.kBT_fullquench = conditions.kB * conditions.T_fullquench;

% logic checks on temperatures
if size(conditions.T_equilibrium,1)~=1
    error('T_equilibrium should be a row vector')
elseif ~isscalar(conditions.T_fullquench)
    error('quench temperature has to be a scalar')
end
%%% end setting temperatures  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% Set Semiconductor Bandstructure and Lattice sites %%%%%%%%%%%%%%%%%%%%%%
% this has to come after the temperature vector is set
%%% Ga2O3 %%

conditions.vibent_T0 = 870; % This is the characteristic To for the average phonon mode in Ga2O3, as determined from the wo that makes the Debeye temperature come out to about 740 K based on the Cv(T) data found online.  To where x = hbar*w0/kBT = T0/T, wherein w0 is the mean phonon energy from the Debeye temperature determined from Cv(T).  

% when T-independent values are used
conditions.TRef = 300;
conditions.EgRef = 4.8;   % this is for constant Eg(T)
conditions.EvRef = 0;
conditions.EcRef = conditions.EgRef;
conditions.NcRef = 1e19;
conditions.NvRef = 5e20;


% T-dependent Eg values
conditions.Eg0 = 5.1;
% conditions.Eg0 = 5.0;
% conditions.Eg0 = 4.9;   % bandgap at 0K in eV.  Using Adrian's number from STEM EELS at PSU - 4.8 eV at 300 K.  So about 4.9 eV at 0K using
conditions.varshini_a = 0.00105248259206626;  % can implement any model you want for Eg(T) as long as it gives a value for all T_equilibrium
conditions.varshini_b = 676.975385507958;

% conditions.EcT_fraction = 0;
% conditions.EcT_fraction = 0.25;
conditions.EcT_fraction = 0.30;
% conditions.EcT_fraction = 0.375;
% conditions.EcT_fraction = 0.50;     % what fraction of Eg(T)  happens in the CB?
% conditions.EcT_fraction = 0.75;
% conditions.EcT_fraction = 1;
conditions.EvT_fraction = 1-conditions.EcT_fraction;

% these take the equilibrium T values from conditions and generate the
% needed T-dependent values.  The fraction of Eg change caused by VB vs CB
% can be modified above.
conditions.NcT_equilibrium = conditions.NcRef*(conditions.T_equilibrium/conditions.TRef).^1.5;
conditions.NvT_equilibrium = conditions.NvRef*(conditions.T_equilibrium/conditions.TRef).^1.5;
delta_EgT_equilibrium = (conditions.varshini_a .* conditions.T_equilibrium.^2)./(conditions.T_equilibrium + conditions.varshini_b);  %delta is a positive number
conditions.EgT_equilibrium = conditions.Eg0 - delta_EgT_equilibrium;
conditions.EcT_equilibrium = conditions.Eg0*ones(size(conditions.EgT_equilibrium)) - conditions.EcT_fraction*delta_EgT_equilibrium;
conditions.EvT_equilibrium = zeros(size(conditions.EgT_equilibrium)) + conditions.EvT_fraction*delta_EgT_equilibrium;

conditions.NcT_fullquench = conditions.NcRef*(conditions.T_fullquench/conditions.TRef).^1.5;
conditions.NvT_fullquench = conditions.NvRef*(conditions.T_fullquench/conditions.TRef).^1.5;
delta_EgT_fullquench = (conditions.varshini_a .* conditions.T_fullquench.^2)./(conditions.T_fullquench + conditions.varshini_b);  %delta is a positive number
conditions.EgT_fullquench = conditions.Eg0 - delta_EgT_fullquench;
conditions.EcT_fullquench = conditions.Eg0 - conditions.EcT_fraction*delta_EgT_fullquench;
conditions.EvT_fullquench = conditions.EvT_fraction*delta_EgT_fullquench;


clear delta_EgT_equilibrium delta_EgT_fullquench   %clean up memory after using these

% site densities in the Ga2O3 lattice (num FU/unit cell / vol for unit cell)
N_Ga1 = 1.91e22;    %site 1
N_Ga2 = 1.91e22;  % site 2
N_O1 = 1.91e22;  %site 3
N_O2 = 1.91e22;   %site 4
N_O3 = 1.91e22;  %site5
N_ia = 1.91e22;
N_ib = 1.91e22;
N_ic = 1.91e22;
N_unspecified_i = 1.91e22;  % interstitials - watch the degeneracy of these
conditions.num_sites = [N_Ga1; N_Ga2; N_O1; N_O2; N_O3; N_ia; N_ib; N_ic; N_unspecified_i];   % note this needs to be a column vector, dont change it to a row vector

% % Calculate the numerical prefactor for each defect from defects.degen_factor_config, defects.degen_factor_elec, defects.cs_num_each_site, and conditions.num_sites

%% defects.cs_num_each_site tells you the number of primitive unit cells needed to form each defect (assuming the sites of the crystal, including distinct interstitials like the ia, ib, ... in b-ga2O3, are counted such that there is one per primitive unit cell
defects.cs_site_prefactor = sort((ones(defects.num_chargestates,1)*conditions.num_sites')./defects.cs_num_each_site,2);      % calculates the max number of each chargestate that could be formed, given the availability of each site type in the crystal (or unit cells, or supercells, whatever the basis for counting).  Then sorts by the rows ascending, so the first column becomes the max number we could form.  Its done with a sort becasue you end up with Inf when you divide by sites you need zero of
defects.cs_site_prefactor = defects.cs_site_prefactor(:,1);  % toss out the other columns keeping only the first one, which should be the limiting one
defects.cs_prefactor = defects.cs_degen_factor_config .* defects.cs_degen_factor_elec .* defects.cs_site_prefactor;

clear N_Ga1 N_Ga2 N_O1 N_O2 N_O3 N_ia N_ib N_ic N_unspecified_i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%% Determine how to set pressures, chemical potentials, concentrations.. of the elements


%% set total pressure and units %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set total pressure
% conditions.P_tot = 1e-4;   % total pressure in atm
conditions.P_tot = 1;   % total pressure in atm
% conditions.P_tot = 0.025;   % total pressure in atm

conditions.P_units = 'atm';  %p_ref set inside each G0 file


% conditions.P_tot = 30e-3;
% conditions.P_units = 'Torr';



%% Set up the elements to use in the model %%%%%%%%%%%%%
% note: later below some elements or individual defects' concentrations can
% be fixed.  This overrides setting chemical potentials up here so it's OK
% to have entries here for elements fixed later.

% conditions.num_elements = 17;
% if conditions.num_elements is not explicitly set, set it to the number of
% columns in defects.cs_dm
if ~isfield(conditions,'num_elements') || sum(isempty(conditions.num_elements))
    conditions.num_elements = size(defects.cs_dm,2);
end


%% Set the scenario used to determine the chemical potential of the host material.
%%% different options for setting the host material's chemical potentials according to equilibrium equations or other

% for T-independent mu conditions (set one mu like mu=0 for rich/poor conditions)
% conditions.mu_Ga2O3 = -10.22;        % formation energy of Ga2O3 in eV/FU from Joel

conditions.mu_conditions_flag = 'pO2-variable';   % equilibrium with pO2 set to arbitrary value.  you must also set pO2.
% conditions.mu_conditions_flag = 'O-rich';
% conditions.mu_conditions_flag = 'Ga-rich';
% conditions.mu_conditions_flag = 'O2-1atm';   % pure O2 gas at 1 atm
% conditions.mu_conditions_flag = 'air-1atm';   % 1 atm of 80/20 N2/O2
% conditions.mu_conditions_flag = 'pGa_vap-variable';   % equilibrium with pGa_vapor set to arbitrary value.  you must also set pGa_vap.
% conditions.mu_conditions_flag = 'Ga_equilib_vap_pressure';    % set pGa_vap to equilibrium vapor pressure of Ga over Ga
% conditions.mu_conditions_flag = 'Ga_ls-rich';  % use this to mean Ga2O3 in direct contact with pure Ga liq or solid
% conditions.mu_conditions_flag = 'pGa2O_vap-variable';
% conditions.mu_conditions_flag = 'pGaO_vap-variable';

% if using "pO2-variable", need to also set pO2
% conditions.pO2 = 1;
% conditions.pO2 = 0.2;
% conditions.pO2 = 0.02;
% conditions.pO2 = 0.2*conditions.P_tot;  % for MOCVD
% conditions.pO2 = 1e-4;
% conditions.pO2 = 1e-6;
% conditions.pO2 = 1e-8;

%% %%%%%%%%%%%%%% Start Looping Here for Pressure %%%%%%%%%%%%%%%%
% Initialize pO2 range and cell array to store outputs of each iteration
conditions.pO2_range = linspace(1e-30,1,10); 
num_pO2 = length(conditions.pO2_range); % Number of pO2 points
% Loop through each value in the range for pO2s

for pIdx = 1:num_pO2
    % Access the current value of pO2 to use in calculations
    conditions.pO2 = conditions.pO2_range(pIdx);

    % Print out current pO2 (debugging)
    fprintf('Current pO2: %e\n', conditions.pO2);


% if using "pGa_vap-variable", need to also set a value
% conditions.pGa_vap = 1e-4;
% conditions.pGa_vap = 1e-8;
% conditions.pGa_vap = 1e-10;

% specify partial pressures of Ga2O or GaO
% conditions.pGa2O_vap = 1e-4;
% conditions.pGaO_vap = 1e-4;

if strcmp(conditions.mu_conditions_flag,'O2-1atm')   % equilibrium with O2 at pO2=ptot=1 atm.  This overrides anything you set accidentally above.
    conditions.pO2 = conditions.pO2(pIdx); % edited for sweeping pressure: we want the partial pressure of oxygen to vary
    conditions.p_tot = 1;
end
%% end setting scenario for host material





%% set chemical potentials depending on the method for each element

% Set default that all elements have fixed value of mu that is
% very negative like -30 eV.  Make a matrix with rows for each T and columns for each
% element (including matrix ones).
conditions.mu_constant_flag = ones(1,conditions.num_elements);

% intialize a flag for mu set by phase boundaries - this should be mutually
% exclusive with the constant mu flag above.
conditions.mu_from_2nd_phases_flag = zeros(1,conditions.num_elements);   % set all to zero/false first then overwrite if needed

% both of these options will overwrite entries in this variable intialized here - this is the one the actual calc uses (so dont change its name unless doing it inside the main calc too)
conditions.muT_equilibrium = -30*ones(conditions.num_T_equilibrium, defects.numelements); %

%% Matrix elements: set mu values first (first 2 columns for a binary, first 3 for a ternary, etc)
[mu_Ga, mu_O] = Ga_O_Ga2O3_chem_potentials_from_conditions(conditions);
conditions.muT_equilibrium(:,1) = mu_Ga;  % these 3 lines not combined to keep it explicit that first 2 columns are Ga and O
conditions.muT_equilibrium(:,2) = mu_O;

%% Impurity elements
% calculate the values of all elements' limiting mu values based on the conditions
% we use only some of them as needed - use them by overwriting the muT_equilibrium
% for Ga2O3, we are using limits based on oxygen, but other limits can also
% be calculated and imposed, sequentially even (O can exclude some ranges for some elements,
% Ga can exclude other ranges for some elements, ..)
[impurity_mu_boundary_limits] = impurity_mu_limit_imposed_by_pO2(conditions.T_equilibrium, conditions.P_tot, conditions.P_units, mu_O);  % note this function's output first column is the first impurity element


% set either constant mu values, T-dependent mu values, or mu values from equilibrium with 2nd phases.  Flip the flag values from 0/1 for fixed mu and mu from 2nd phases

% example for setting Si to -4 for all T.  the const mu flag is already 1 and from_2nd_phases flag =0 since this is how they're intialized but explicitly set here for clarity
% conditions.mu_constant_flag(3) = 1;
% conditions.mu_from_2nd_phases_flag(3) = 0
% mu_Si = -4;
% conditions.muT_equilibrium(:,3) = mu_Si*ones(conditions.num_T_equilibrium,1);

% % example for setting H to value set by 2nd phase equilibrium based on mu_O
% H_mu_num = 4;
% conditions.mu_from_2nd_phases_flag(H_mu_num) = 1;
% conditions.mu_constant_flag(H_mu_num) = 0;   % flip the other flag too
% conditions.muT_equilibrium(:,H_mu_num) = impurity_mu_boundary_limits(:,H_mu_num-2);  % Use the right precomputed column.  Column numbers for impurity_mu variable are smaller by 2 because first 2 columns are Ga and O

% % mu_Sn set by 2nd phase equilibrium based on mu_O
% Sn_mu_num = 6;
% conditions.mu_from_2nd_phases_flag(Sn_mu_num) = 1;
% conditions.mu_constant_flag(Sn_mu_num) = 0;   % flip the other flag too
% conditions.muT_equilibrium(:,Sn_mu_num) = impurity_mu_boundary_limits(:,Sn_mu_num-2);  % Use the right precomputed column.  Column numbers for impurity_mu variable are smaller by 2 because first 2 columns are Ga and O


% logic checks
mu_set_exclusively_by_one_method_check = ((conditions.mu_constant_flag + conditions.mu_from_2nd_phases_flag)-1) ~=0 ;
conditions.mu_set_flag = (conditions.mu_constant_flag + conditions.mu_from_2nd_phases_flag)>0;

if any(mu_set_exclusively_by_one_method_check)
    disp(mu_set_exclusively_by_one_method_check)
    error('mu for each element has to be set either as a constant or by a phase boundary.  The line above has a 1 in the position of the element causing the problem')
end

%%%%%%% finished specifying all mu values for elements controlled by mu rather than total conc.



%% Set elements with fixed concentrations and the target values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if an element is fixed, it means the chemical potential for that element
% is floated (rather than fixed) and optimized (along with Ef and those of
% any others fixed )to get the desired total number over all the defects
% containing that element.  List any fixed elements by putting a 1 in its
% position in the vector conditions.fixed_conc_flag.
% This is mutually exclusive with setting chem potential for the element
%
% element order for ref
% 1=Ga 2=O 3=Si 4=H 5=Fe 6=Sn 7=Cr 8=Ti 9=Ir 10=Mg 11=Ca 12=Zn 13=Co 14=Zr
% 15=Hf 16=Ta 17=Ge, 18=Pt 19= Rh


% initilaize flag and value pair for elements set by total concentrations -
% since this is a more unusual option this supercedes any mu values set (mu values determined in the calc self consistently).  This should be mutually exclusive with setting mu value flags though.
conditions.fixed_conc_flag = zeros(1,conditions.num_elements);
conditions.fixed_conc_values = zeros(1,conditions.num_elements);  % this will hold target values of elements

% initialize holder for the ranges over which to search for mu for frozen elements.  convention for now is that this matrix should be the size for all the elements, not just the fixed ones we specify ranges for all the elements.  Later can change this to only specify ranges for the ones that are frozen (not specify values for those left open).
% since either a direct grid search or particle swarm is used to generate the first guess, there are
% HUGE time savings if you can narrow down this range manually first using
% a series of calcs with constant mu values to estimate mu for each element
% you want to freeze.  Or you can brute force it and just wait for ever...%
conditions.fixed_elements_mu_ranges = zeros(conditions.num_elements,2);   % create  matrix of right size for ALL the elements - only will pay attention to those which are fixed


% 
% conditions.fixed_conc_flag(3) = 1;          % fix Si concentration
% conditions.mu_set_flag(3) = 0;              % flip the flag saying Si is set by mu
% conditions.fixed_conc_values(3) = 5e16;     % specify the value
% lo_mu_Si = -10; % set range to search for the unknown mu
% hi_mu_Si = -2;
% conditions.fixed_elements_mu_ranges(3,:) = [lo_mu_Si hi_mu_Si];
% 

%
% %% loop over si doping -
%
% Si_doping = [1e15 1e16 1e17 1e18 1e19];
% conditions.fixed_conc_flag(3) = 1;          % fix Si concentration
% conditions.mu_set_flag(3) = 0;              % flip the flag saying Si is set by mu
%
%
% for index = 1:numel(Si_doping)
%
%     save_fname_root = strcat('Si_OMVPE_',num2str(Si_doping(index)));
%     conditions_save_fname = strcat(save_fname_root,'_calculation_conditions');
%     equilib_save_fname = strcat(save_fname_root,'_equilib_dark_sol');
%     fullquench_save_fname = strcat(save_fname_root,'_fullquench_dark_sol');
%
%     conditions.fixed_conc_values(3) = Si_doping(index);
%

% conditions.fixed_conc_flag(3) = 1;    % fix Si
% conditions.mu_set_flag(3) = 0;
% conditions.fixed_conc_values(3) = 1e17;
% lo_mu_Si = -9;
% hi_mu_Si = -6;
% conditions.fixed_elements_mu_ranges(3,:) = [lo_mu_Si hi_mu_Si];


% conditions.fixed_conc_flag(4) = 1;    % fix H
% conditions.mu_set_flag(4) = 0;
% conditions.fixed_conc_values(4) = 1e17;
% % H is about -2 to -4 eV for 4.5e18 Sn and 1e17 H
% lo_mu_H = -4;
% hi_mu_H = -2;
% conditions.fixed_elements_mu_ranges(4,:) = [lo_mu_H hi_mu_H];

% % conditions.fixed_conc_flag(5) = 1;    % fix Fe
% % conditions.mu_set_flag(5) = 0;
% % conditions.fixed_conc_values(5) = 2.8e16;
% % % % 1e16 to 1e17 Fe is about -6 to -5 eV with 2e18 Sn doping
% % lo_mu_Fe = -7;
% % hi_mu_Fe = -4;
% % conditions.fixed_elements_mu_ranges(5,:) = [lo_mu_Fe hi_mu_Fe];
% 
% conditions.fixed_conc_flag(6) = 1;    % fix Sn
% conditions.mu_set_flag(6) = 0;              % flip the flag saying Sn is set by mu
% conditions.fixed_conc_values(6) = 4.5e18;
% % % 2e18 Sn is about -5.5 eV at 2100K  to -5 eV at 500 K
% lo_mu_Sn = -6.5; % set range to search for the unknown mu
% hi_mu_Sn = -3;
% conditions.fixed_elements_mu_ranges(6,:) = [lo_mu_Sn hi_mu_Sn];

% conditions.fixed_conc_flag(7) = 1;    % fix Cr
% conditions.mu_set_flag(7) = 0;
% conditions.fixed_conc_values(7) = 1.9e16;
% % 1e16 Cr is about -8 to -9 for 500-2100 K
% lo_mu_Cr = -10;
% hi_mu_Cr = -4;
% conditions.fixed_elements_mu_ranges(7,:) = [lo_mu_Cr hi_mu_Cr];

% conditions.fixed_conc_flag(8) = 1;    % fix Ti at 1e16
% conditions.mu_set_flag(8) = 0;
% conditions.fixed_conc_values(8) = 1e16;
% % 1e16 Ti
% lo_mu_Ti = -10;
% hi_mu_Ti = -5;
% conditions.fixed_elements_mu_ranges(8,:) = [lo_mu_Ti hi_mu_Ti];

% conditions.fixed_conc_flag(9) = 1;    % fix Ir at 1.7e17
% conditions.mu_set_flag(9) = 0;
% conditions.fixed_conc_values(9) = 1.7e17;
% lo_mu_Ir = -3;
% hi_mu_Ir = -1.5;
% conditions.fixed_elements_mu_ranges(9,:) = [lo_mu_Ir hi_mu_Ir];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% the if loop below implements a batch loop to make 8 scenarios, permuting Fe, Sn, Mg,
% and UID doping with f=0.5, 0.75.  "index" is the index variable used in
% overall for loop
%
% if index==1 || index==5
%
%     % Sn KURAMATA 2016
    % conditions.fixed_conc_flag(3) = 1;    % fix Si
    % conditions.mu_set_flag(3) = 0;              % flip the flag saying Si is set by mu
    % conditions.fixed_conc_values(3) = 2.0e17;
    % conditions.fixed_conc_flag(6) = 1;    % fix Sn
    % conditions.mu_set_flag(6) = 0;
    % conditions.fixed_conc_values(6) = 4.5e18;
    % conditions.fixed_conc_flag(5) = 1;    % fix Fe
    % conditions.mu_set_flag(5) = 0;
    % conditions.fixed_conc_values(5) = 6.4e16;
%
% elseif index==2 || index==6
%
    % % Mg
    % conditions.fixed_conc_flag(3) = 1;    % fix Si
    % conditions.mu_set_flag(3) = 0;
    % conditions.fixed_conc_values(3) = 2.3e17;
    % % conditions.fixed_conc_flag(6) = 1;    % fix Sn
    % % conditions.mu_set_flag(6) = 0;
    % % conditions.fixed_conc_values(6) = 6.3e15;
    % % conditions.fixed_conc_flag(5) = 1;    % fix Fe
    % % conditions.mu_set_flag(5) = 0;
    % % conditions.fixed_conc_values(5) = 3.0e16;

% elseif index==3 || index==7
%
% Fe
% conditions.fixed_conc_flag(3) = 1;    % fix Si
%     conditions.mu_set_flag(3) = 0;
% conditions.fixed_conc_values(3) = 2.3e17;
% conditions.fixed_conc_flag(6) = 1;    % fix Sn
%     conditions.mu_set_flag(6) = 0;
% conditions.fixed_conc_values(6) = 6.5e15;
% conditions.fixed_conc_flag(5) = 1;    % fix Fe
%     conditions.mu_set_flag(5) = 0;
% conditions.fixed_conc_values(5) = 2.5e18;
%
% elseif index==4 || index==8
%     % UID
%     conditions.fixed_conc_flag(3) = 1;    % fix Si
%     conditions.mu_set_flag(3) = 0;
%     conditions.fixed_conc_values(3) = 2.3e17;
%     conditions.fixed_conc_flag(6) = 1;    % fix Sn
%     conditions.mu_set_flag(6) = 0;
%     conditions.fixed_conc_values(6) = 6.5e15;
%     conditions.fixed_conc_flag(5) = 1;    % fix Fe
%     conditions.mu_set_flag(5) = 0;
%     conditions.fixed_conc_values(5) = 3.0e16;
%
% end


% to do: add in ability to import T-dependnet Ef and mu's from a prior simulation.  Let's say you fixed 2 elements, got a good solution, and next want to freeze a 3rd one.  The best place to start is that last solution then ramp up the fixed concentration on the new element from 0 to it's target value.
% actually could make it automatically do this for multi-element problems - start with the largest target concentration element, then add the next a perturbing one especially if it forms complexes with that first one, then go to the next element..

% default is no elements fixed - catch if nothing entred.
if ~isfield(conditions,'fixed_conc_flag') || sum(isempty(conditions.fixed_conc_flag))
    conditions.fixed_conc_flag = zeros(defects.numelements,1);
elseif ~isfield(conditions,'fixed_conc_values') || sum(isempty(conditions.fixed_conc_values))
    conditions.fixed_conc_values = zeros(defects.numelements,1);
end

% if something is not fixed but by accident it has a target conc set, jsut
% set its target concentration to zero also, using the fixed elements vector
% as a mask
conditions.fixed_conc_values = conditions.fixed_conc_values.*conditions.fixed_conc_flag;

% find number of fixed elements and their indices
conditions.indices_of_fixed_elements = find(conditions.fixed_conc_flag);
conditions.num_fixed_elements = numel(conditions.indices_of_fixed_elements);

if conditions.num_fixed_elements~=0
    for i=1:conditions.num_fixed_elements
        msg = strcat('WARNING: Fixed concentration of element used for...',defects.elementnames(conditions.indices_of_fixed_elements(i)),'...Be sure you meant to do this.');
        disp(msg)
    end
elseif conditions.num_fixed_elements==0
    disp('No Eelements have fixed concentrations.  Proceeding')
end

clear i msg impurity_mu_boundary_limits mu_set_exclusively_by_one_method_check max_mu_range start_kBT fine_factor
clear mu_Ga mu_O lo_mu_Ga hi_mu_Ga lo_mu_O hi_mu_O lo_mu_Si hi_mu_Si lo_mu_H hi_mu_H lo_mu_Fe hi_mu_Fe lo_mu_Sn hi_mu_Sn lo_mu_Cr hi_mu_Cr lo_mu_Ti hi_mu_Ti lo_mu_Ir hi_mu_Ir
%%% end section on fixed elements %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% logic checks on mu or concentraitons set for each element
% catch it if by accident the mu and defect conc are both set for any
% elements

neither_mu_nor_conc_flag = (conditions.mu_set_flag .* conditions.fixed_conc_flag);
both_mu_and_conc_flag = conditions.mu_set_flag + conditions.fixed_conc_flag ==0;

if any(neither_mu_nor_conc_flag)
    disp(neither_mu_nor_conc_flag)
    error('Must set either the mu or concentation for each element.')
elseif any(both_mu_and_conc_flag)
    disp(both_mu_and_conc_flag)
    error('Must set either the mu or the concentation for each element, not both.')
elseif size(conditions.muT_equilibrium,2)~=conditions.num_elements
    error('number of columns in mu array not equal to number of elements declared')
elseif sum(sum(isinf(conditions.muT_equilibrium)))
    error('Infinite mu detected - this happens when there is no model for G0 for the T requested for some substance(s)')
end

clear both_mu_and_conc_flag neither_mu_nor_conc_flag
%%% end handling matrix chemical potentials %%%%%%%%%%%%




%% set unidentified shallow doping (fixed charge added to charge balance)
%conditions.Nd = 0;   % imposed shallow doping concentrations
%conditions.Nd = 2e16;   % UID in real samples with Si backgroud
% conditions.Nd = 2e17;   % intentional doping
% conditions.Nd = 0;
% conditions.Na = 0;

% set shallow dopants to zero as default
if ~isfield(conditions,'Nd') || isempty(conditions.Nd)
    conditions.Nd = 0;
end

if ~isfield(conditions,'Na') || isempty(conditions.Na)
    conditions.Na = 0;
end
%%% end shallow doping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%
% %% set any interstitials that equilibrate at T_quench
% if strcmp(conditions.interstitial_equilibrate_at_Tquench_flag, 'On')
%     cs_is_simple_interstitial_flag = ~any(defects.cs_num_each_site(:,1:5)) && sum(defects.cs_dm,2)==1;  % if the cs uses none of the {Ga1 Ga2 O1 O2 O3} sites, thus uses only interstital sites AND it adds only one atom, it's a simple interstitial.
%     % also add other defects or chargstates explicitly by manually calling
%     % out their index number here, or a different logic condition
%
% elseif strcmp(conditions.interstitial_equilibrate_at_Tquench_flag, 'Off')
%     cs_is_simple_interstitial_flag = zeros(defects.num_chargestates,1);
% else
%     error('Something weird about some chargestates in terms of whether or not they are simple interstitials')
% end
%
%
% for kk = 1:defects.num_chargestates
%     if  cs_is_simple_interstitial_flag(kk)
%         conditions.set_to_zero_for_paraequilibrium_flag = 1;
%     % elseif some other condition based on diffusivity or, it includes small element like Li or H ...
%     end
% end







%% Set any defects with fixed concentrations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this overrides any conditions set for the elements.
% first initialize a COLUMNvector with the right number of entries for the
% current defect variable.  Then you can go in and st individual ones
conditions.fixed_defects = zeros(defects.num_defects,1);
conditions.fixed_defects_concentrations = zeros(defects.num_defects,1);

% %% zero out the   % VGa-2Sn and -3Sn
% conditions.fixed_defects_concentrations(163:169) = 0; 
% conditions.fixed_defects(163:166)= 1;






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% T independent case - each defect has single value for all temperatures

%
% conditions.fixed_defects(169)= 1;
% % conditions.fixed_defects_concentrations(169) = 8.9e16;   % Si_GaI
% %
% %  %% Sn doped Kuramata 2016
%     conditions.fixed_defects(222)= 1;
%     conditions.fixed_defects_concentrations(222) = 3.6e17;   % Ir_GaII
%     conditions.fixed_defects(193)= 1;
%     conditions.fixed_defects_concentrations(193) = 4.3e16;   % Cr_GaII
%     conditions.fixed_defects(199)= 1;
%     conditions.fixed_defects_concentrations(199) = 2.5e15;   % Ti_GaII
%     conditions.fixed_defects(247)= 1;
%     conditions.fixed_defects_concentrations(247) = 6.3e15;   % Zr_GaII
%     conditions.fixed_defects(238)= 1;
%     conditions.fixed_defects_concentrations(238) = 1.3e16;   % Ca_GaII
%     conditions.fixed_defects(236)= 1;
%     conditions.fixed_defects_concentrations(236) = 1.4e16;   % Mg_GaII
% %
% if index==1 || index==5
%     %% Sn doped Kuramata 2016
%     conditions.fixed_defects(222)= 1;
%     conditions.fixed_defects_concentrations(222) = 3.6e17;   % Ir_GaII
%     conditions.fixed_defects(193)= 1;
%     conditions.fixed_defects_concentrations(193) = 4.3e16;   % Cr_GaII
%     conditions.fixed_defects(199)= 1;
%     conditions.fixed_defects_concentrations(199) = 2.5e15;   % Ti_GaII
%     conditions.fixed_defects(247)= 1;
%     conditions.fixed_defects_concentrations(247) = 6.3e15;   % Zr_GaII
%     conditions.fixed_defects(238)= 1;
%     conditions.fixed_defects_concentrations(238) = 1.3e16;   % Ca_GaII
%     conditions.fixed_defects(236)= 1;
%     conditions.fixed_defects_concentrations(236) = 1.4e16;   % Mg_GaII
%
% elseif index==2 || index==6
% %
%     %% Mg doped
%     conditions.fixed_defects(222)= 1;
%     conditions.fixed_defects_concentrations(222) = 5.9e16;   % Ir_GaII
%     conditions.fixed_defects(193)= 1;
%     conditions.fixed_defects_concentrations(193) = 3.9e16;   % Cr_GaII
%     conditions.fixed_defects(199)= 1;
%     conditions.fixed_defects_concentrations(199) = 3.8e15;   % Ti_GaII
%     conditions.fixed_defects(247)= 1;
%     conditions.fixed_defects_concentrations(247) = 8.8e15;   % Zr_GaII
%     conditions.fixed_defects(238)= 1;
%     conditions.fixed_defects_concentrations(238) = 1.3e16;   % Ca_GaII
%     conditions.fixed_defects(236)= 1;
%     conditions.fixed_defects_concentrations(236) = 2e18;   % Mg_GaII
% %
% elseif index==3  || index==7
%     %% Fe doped
%     conditions.fixed_defects(222)= 1;
%     conditions.fixed_defects_concentrations(222) = 5.9e16;   % Ir_GaII
%     conditions.fixed_defects(193)= 1;
%     conditions.fixed_defects_concentrations(193) = 3.9e16;   % Cr_GaII
%     conditions.fixed_defects(199)= 1;
%     conditions.fixed_defects_concentrations(199) = 3.8e15;   % Ti_GaII
%     conditions.fixed_defects(247)= 1;
%     conditions.fixed_defects_concentrations(247) = 8.8e15;   % Zr_GaII
%     conditions.fixed_defects(238)= 1;
%     conditions.fixed_defects_concentrations(238) = 1.3e16;   % Ca_GaII
%     conditions.fixed_defects(236)= 1;
%     conditions.fixed_defects_concentrations(236) = 5e15;   % Mg_GaII

% elseif index==4  || index==8
%
%     %% UID Kuramata 2016
%     conditions.fixed_defects(222)= 1;
%     conditions.fixed_defects_concentrations(222) = 5.9e16;   % Ir_GaII
%     conditions.fixed_defects(193)= 1;
%     conditions.fixed_defects_concentrations(193) = 3.9e16;   % Cr_GaII
%     conditions.fixed_defects(199)= 1;
%     conditions.fixed_defects_concentrations(199) = 3.8e15;   % Ti_GaII
%     conditions.fixed_defects(247)= 1;
%     conditions.fixed_defects_concentrations(247) = 8.8e15;   % Zr_GaII
%     conditions.fixed_defects(238)= 1;
%     conditions.fixed_defects_concentrations(238) = 1.3e16;   % Ca_GaII
%     conditions.fixed_defects(236)= 1;
%     conditions.fixed_defects_concentrations(236) = 5e15;   % Mg_GaII
%
% end
%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T dependent fixed defect values - use this for example after floating the
% element and finding that only 2 defects dominate but their numbers change
% vs T.  The format will be T down the rows and concentrations for each
% defect across.  So to pick out the column vector giving concentration of
% defect #32 at all the temperatures you'd use conditions.fixed_defects_concentrations(:,32)
%
if strcmp(conditions.T_dep_fixed_defect_flag,'On') && size(conditions.fixed_defects_concentrations,2)==1   % if the T dependent flag is on and the list of fixed defects so far is just a column vector with single values
    conditions.fixed_defects_concentrations = ones(size(conditions.T_equilibrium,2),1) * conditions.fixed_defects_concentrations';  % expand out the list so each T can have a distinct value.

    % now we can go in and assign different values for different defects, for
    % example based on values obtained from a prior simuation.  Yes for now specify the fixed defects
    % inside the if loop checking if the flag is on

    %     conditions.fixed_defects(15)= ones(1,1);  % VGa,ic in the 3/21/24 database
    %% this one for 2100:-20:500
    % conditions.fixed_defects_concentrations(:,15) = [2.26508E+11; 2.01586E+11; 1.85099E+11; 1.76221E+11; 1.74721E+11; 1.80022E+11; 1.96051E+11; 2.218E+11; 2.61093E+11; 3.18136E+11; 3.96893E+11; 5.12249E+11; 6.7033E+11; 8.90973E+11; 1.1998E+12; 1.634E+12; 2.24789E+12; 3.12142E+12; 4.37307E+12; 6.17981E+12; 8.80822E+12; 1.26631E+13; 1.83646E+13; 2.68714E+13; 3.96789E+13; 5.91426E+13; 8.90082E+13; 1.35292E+14; 2.07752E+14; 3.22372E+14; 5.05573E+14; 8.01379E+14; 1.28E+15; 2.08E+15; 3.38E+15; 5.54E+15; 9.03E+15; 1.45E+16; 2.27E+16; 3.37E+16; 4.66E+16; 6.12E+16; 7.31E+16; 8.27E+16; 8.85E+16; 9.05E+16; 8.88E+16; 8.41E+16; 7.71E+16; 6.87E+16; 5.98E+16; 5.03E+16; 4.14E+16; 3.31E+16; 2.60E+16; 1.98E+16; 1.46E+16; 1.07E+16; 7.60E+15; 5.11E+15; 3.43E+15; 2.20E+15; 1.36E+15; 8.11042E+14; 4.64816E+14; 2.55007E+14; 1.34026E+14; 6.59593E+13; 3.22555E+13; 1.47954E+13; 6.45521E+12; 2.69887E+12; 8.25486E+11; 3.0156E+11; 1.10664E+11; 34393838843; 10298916679; 2931857896; 768916197.3; 183911746.7; 39640811.58];
    %% this one for 1500:-100:500
    %     conditions.fixed_defects_concentrations(:,15) = [5.05573E+14; 5.54E+15; 4.66E+16; 9.05E+16; 5.98E+16; 1.98E+16; 3.43E+15; 2.55007E+14; 6.45521E+12];

end   %end the if loop checking if T dependent fixed defects are enabled











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% logic checks on fixed defects and setting default to none fixed
% default is none, catch this.  Note that if elseif loops exit once one
% condition is satisfied

if ~isfield(conditions,'fixed_defects') || sum(isempty(conditions.fixed_defects))
    conditions.fixed_defects = zeros(defects.num_defects,1);
end

if strcmp(conditions.T_dep_fixed_defect_flag,'Off')
    conditions.fixed_defects_concentrations = conditions.fixed_defects_concentrations.*conditions.fixed_defects;          % just in case a defect is not fixed but by accident it has a target conc set, just set its target concentration to zero also, using the fixed defects vector as a mask

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
    for i=1:conditions.num_fixed_defects
        msg = strcat('WARNING: Fixed concentration used for...',defects.defect_names(conditions.indices_of_fixed_defects(i)),'...Be sure you meant to do this.');
        disp(msg)
    end
elseif conditions.num_fixed_defects==0
    disp('No Defects fixed.  Proceeding')
end

%%%  end section on fixed defects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%% save the problem definition in the conditions variable into a mat file for later inspection
save(strcat(conditions.save_pname,conditions.conditions_save_fname,'.mat'),'conditions')
disp('Saved problem specificaiton conditions to .mat file')
% clear conditions_save_fname
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% end problem definition and calc details section %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% some automatic input validation before proceeding %%%%%%%%%%%%%%%%%%%%%%
if size(conditions.muT_equilibrium(1,:),2) ~= size(defects.cs_dm,2)
    error('number of rows in mu must equal number of columns in cs_dm')
elseif (size(defects.cs_dHo,1) ~= size(defects.cs_charge,1))
    error('dE and cs_charge must be the same size')
    % elseif size(conditions.num_sites,1) ~= size(defects.cs_dm,2)
    %     error('sites is a column vector giving the number of lattice sites.  cs_dm is a matrix with rows for each defect and number of atoms added or subtracted from each site in the columns')
elseif max(size(conditions.T_fullquench))~=1
    error('for now, T_quench needs to be a scalar - vectorize it later?')
    % elseif size(defects.gq_isfrozen_allT_defects,1)~= defects.num_defects || size(defects.gq_isfrozen_allT_concdefects,1)~= defects.num_defects || size(defects.gq_isfrozen_allT_chargestates,1)~=defects.num_chargestates
    %     error('check sizes of defects.gq_isfrozen_allT_defects, defects.gq_isfrozen_allT_concdefects, and defects.gq_isfrozen_allT_chargestates - should be column vectors same size as number of defects')
    % elseif any((defects.gq_isfrozen_allT_defects==1) ~= (defects.gq_isfrozen_allT_concdefects~=0))
    %     error('if a defect is frozen for all temperatures by placing a 1 in defects.gq_isfrozen_allT_defects, then a concentration mut be given for it in defects.gq_isfrozen_allT_concdefects')
    % elseif size(defects.gq_isfrozen_allT_chargestates_index,1)+size(defects.gq_isNOTfrozen_allT_chargestates_index,1)~=defects.num_chargestates
    %     error('something wrong with numbers of frozen and unfrozen chargestates')
    % elseif size(defects.gq_isfrozen_allT_defects_index,1)+size(defects.gq_isNOTfrozen_allT_defects_index,1)~=defects.num_defects
    %     error('something wrong with numbers of frozen and unfrozen defects')
else
    disp('Inputs passed some basic checks on compatibility.  User MUST still check that they did what they really anted to.  Garbage In = Garbagee Out!')
end



%
% %% Plot Formation enthalpies vs EF at 300 K
% figure(1)
% EF = 0:conditions.EgRef/100:conditions.EgRef;
% dH_EF0 = defects.cs_dHo - defects.cs_dm*conditions.muT_equilibrium(1,:)' ;
% for i=1:defects.num_chargestates  % plotting the Ga vacancies only
%     f1line(i) = plot(EF, (dH_EF0(i) + defects.cs_charge(i)*EF),'DisplayName', defects.chargestate_names(i));
%     f1line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f1line(i).DisplayName},size(f1line(i).XData)));
% end
% % for i=1:4  % plotting thedH for Si donors
% % 	plot(EF, (dH_EF0(i) + defects.cs_charge(i)*EF),'DisplayName', defects.chargestate_names(i));
% % end
%
% legend show
% title(strcat('300 K Formation Enthalpies of V_{Ga} for \mu_{O}=',num2str(conditions.muT_equilibrium(1,2)),' and \mu_{Ga}=',num2str(conditions.muT_equilibrium(1,1))))
% % title(strcat('Formation Enthalpies of Si_{Ga} for \mu_{O}=',num2str(mu_O),' and \mu_{Ga}=',num2str(mu_Ga)))
% xlabel('EF (eV)')
% ylabel('Formation Enthalpy (eV)')
% datacursormode.Enable='on';   %turn on the datatip to give name of defect on mouse hover
%
% clear EF dH_EF0 i datacursormode f1line
% %%%%%%%%%%% end formation enthalpy plotting
%




%% Task Trange 1: compute defect equilibrium for a vector of T values
tic
[equilib_dark_sol] = defect_equilibrium_dark(conditions, defects);
toc



% element concentrations.  These are same for equilib and quenching
% element_totals_headers = ["T_equilibrium (K)" "Ga (#/cm3)" "O (#/cm3)" "Si (#/cm3)" "H (#/cm3)" "Fe (#/cm3)" "Sn (#/cm3)" "Cr (#/cm3)" "Ti (#/cm3)" "Ir (#/cm3)" "Mg (#/cm3)" "Ca a(#/cm3)" "Zn (#/cm3)" "Co (#/cm3)" "Zr (#/cm3)" "Hf (#/cm3)" "Ta (#/cm3)" "Ge (#/cm3)"];  % when it is all text entered by hand with no variables concatenated, have to use " " not ' '
element_totals_headers = ['T_equilibrium (K)' defects.elementnames];
element_totals_cell = [cellstr(element_totals_headers); num2cell([equilib_dark_sol.T_equilibrium equilib_dark_sol.element_totals])];
writecell(element_totals_cell,strcat(conditions.save_pname,conditions.elements_save_fname),'Delimiter','tab','WriteMode','overwrite');



%%%%%%%%%%% save the outputs
all_defects_out = cat(2,equilib_dark_sol.T_equilibrium,conditions.pO2,conditions.P_tot,equilib_dark_sol.charge_bal_err, equilib_dark_sol.element_bal_err, equilib_dark_sol.tot_bal_err, equilib_dark_sol.mu, equilib_dark_sol.element_totals, conditions.EvT_equilibrium', conditions.EcT_equilibrium', equilib_dark_sol.EFn, equilib_dark_sol.EFp, equilib_dark_sol.Nd, equilib_dark_sol.Na, equilib_dark_sol.n, equilib_dark_sol.p, equilib_dark_sol.sth1, equilib_dark_sol.sth2, equilib_dark_sol.defects);
all_chargestates_out = cat(2,equilib_dark_sol.T_equilibrium,conditions.pO2,conditions.P_tot,equilib_dark_sol.charge_bal_err, equilib_dark_sol.element_bal_err, equilib_dark_sol.tot_bal_err, equilib_dark_sol.mu, equilib_dark_sol.element_totals, conditions.EvT_equilibrium', conditions.EcT_equilibrium', equilib_dark_sol.EFn, equilib_dark_sol.EFp, equilib_dark_sol.Nd, equilib_dark_sol.Na, equilib_dark_sol.n, equilib_dark_sol.p, equilib_dark_sol.sth1, equilib_dark_sol.sth2, equilib_dark_sol.chargestates);
all_defect_headers = ['T_equilibrium (K)' 'Partial Pressure (atm)' 'Total Pressure (atm)' 'charge_bal_err' 'element_bal_err' 'tot_bal_err' defects.mu_names_with_units defects.elementnames 'Ev(T) (eV)' 'Ec(T) (eV)' 'Efn(eV)' 'Efp (eV)' 'Nd (#/cm3)' 'Na (#/cm3)' 'n (#/cm3)' 'p (#/cm3)' 'sth1 (#/cm3)' 'sth2 (#/cm3)' defects.defect_names'];
all_chargestate_headers = ['T_equilibrium (K)' 'Partial Pressure (atm)' 'Total Pressure (atm)' 'charge_bal_err' 'element_bal_err' 'tot_bal_err' defects.mu_names_with_units defects.elementnames 'Ev(T) (eV)' 'Ec(T) (eV)' 'Efn(eV)' 'Efp (eV)' 'Nd (#/cm3)' 'Na (#/cm3)' 'n (#/cm3)' 'p (#/cm3)' 'sth1 (#/cm3)' 'sth2 (#/cm3)' defects.chargestate_names'];
all_defects_cell = [cellstr(all_defect_headers); num2cell(all_defects_out)];
all_chargestates_cell = [cellstr(all_chargestate_headers); num2cell(all_chargestates_out)];
all_elements_out = cat(2,equilib_dark_sol.T_equilibrium,conditions.pO2,conditions.P_tot,equilib_dark_sol.element_totals);

% throw out all the ones with insignificant concentrations and save
sig_defects_index = find(max(equilib_dark_sol.defects,[],1)>=conditions.save_min);
sig_chargestates_index = find(max(equilib_dark_sol.chargestates,[],1)>=conditions.save_min);
sig_defects_out = cat(2,equilib_dark_sol.T_equilibrium,conditions.pO2,conditions.P_tot, equilib_dark_sol.charge_bal_err, equilib_dark_sol.element_bal_err, equilib_dark_sol.tot_bal_err, equilib_dark_sol.mu, conditions.EvT_equilibrium', conditions.EcT_equilibrium', equilib_dark_sol.EFn, equilib_dark_sol.EFp, equilib_dark_sol.Nd, equilib_dark_sol.Na, equilib_dark_sol.n, equilib_dark_sol.p, equilib_dark_sol.sth1, equilib_dark_sol.sth2, equilib_dark_sol.defects(:,sig_defects_index));
sig_chargestates_out = cat(2,equilib_dark_sol.T_equilibrium,conditions.pO2,conditions.P_tot, equilib_dark_sol.charge_bal_err, equilib_dark_sol.element_bal_err, equilib_dark_sol.tot_bal_err, equilib_dark_sol.mu, conditions.EvT_equilibrium', conditions.EcT_equilibrium', equilib_dark_sol.EFn, equilib_dark_sol.EFp, equilib_dark_sol.Nd, equilib_dark_sol.Na, equilib_dark_sol.n, equilib_dark_sol.p, equilib_dark_sol.sth1, equilib_dark_sol.sth2, equilib_dark_sol.chargestates(:,sig_chargestates_index));
sig_defect_headers = ['T_equilibrium (K)' 'Partial Pressure (atm)' 'Total Pressure (atm)' 'charge_bal_err' 'element_bal_err' 'tot_bal_err' defects.mu_names_with_units 'Ev(T) (eV)' 'Ec(T) (eV)' 'Efn(eV)' 'Efp (eV)' 'Nd (#/cm3)' 'Na (#/cm3)' 'n (#/cm3)' 'p (#/cm3)' 'sth1 (#/cm3)' 'sth2 (#/cm3)' defects.defect_names(sig_defects_index)'];
sig_chargestate_headers = ['T_equilibrium (K)' 'Partial Pressure (atm)' 'Total Pressure (atm)' 'charge_bal_err' 'element_bal_err' 'tot_bal_err' defects.mu_names_with_units 'Ev(T) (eV)' 'Ec(T) (eV)' 'Efn(eV)' 'Efp (eV)' 'Nd (#/cm3)' 'Na (#/cm3)' 'n (#/cm3)' 'p (#/cm3)' 'sth1 (#/cm3)' 'sth2 (#/cm3)' defects.chargestate_names(sig_chargestates_index)'];
sig_defects_cell = [cellstr(sig_defect_headers); num2cell(sig_defects_out)];
sig_chargestates_cell = [cellstr(sig_chargestate_headers); num2cell(sig_chargestates_out)];

% write output tables to a tab delim text file
% writematrix(equilib_dark_sol_out,strcat(conditions.save_pname,conditions.equilib_save_fname),'Delimiter','tab','WriteMode','overwrite')
if strcmp(conditions.save_files_flag,'All')
    writecell(all_defects_cell,strcat(conditions.save_pname,conditions.equilib_save_fname,'_all_defects'),'Delimiter','tab','WriteMode','overwrite');  % write results with column headers
    writecell(all_chargestates_cell,strcat(conditions.save_pname,conditions.equilib_save_fname,'_all_chargestates'),'Delimiter','tab','WriteMode','overwrite');
elseif strcmp(conditions.save_files_flag,'Sig_Only')
    writecell(sig_defects_cell,strcat(conditions.save_pname,conditions.equilib_save_fname,'_sig_defects'),'Delimiter','tab','WriteMode','overwrite');
    writecell(sig_chargestates_cell,strcat(conditions.save_pname,conditions.equilib_save_fname,'_sig_chargestates'),'Delimiter','tab','WriteMode','overwrite');
end
disp('Wrote equilibrium dark solution as four text files')

% write the whole solution structure to a mat file
save(strcat(conditions.save_pname,conditions.equilib_save_fname))
disp('Wrote equilibrium dark solution to .mat file')


%% Storing the Pressure Sweeps for All Defects
    
% Creating a cells to store pressure sweep in the first time through
%   the loop: 
if pIdx == 1
    % Pre-allocate the cell
    num_rows = length(conditions.pO2_range) + 1; % Include the header row
    num_columns = length(all_defect_headers);   % Match number of headers
    all_defect_pO2_sweep_cell = cell(num_rows, num_columns); % Pre-allocate correctly
    % Populate the first row with headers
    all_defect_pO2_sweep_cell(1, :) = cellstr(all_defect_headers); % Store headers
end

% Populate the rest of the rows over each iteration
all_defect_pO2_sweep_cell(pIdx + 1, :) = num2cell(all_defects_out); % Assign data row-wise

%% Storing the Pressure Sweeps for Significant Defects
% 
% % Creating a cell to store each iteration of pressure sweep
% if pIdx == 1
%     % Pre-allocate the cell
%     num_rows = length(conditions.pO2_range) + 1; % Include the header row
%     num_columns = length(sig_defects_out);   % Match number of headers
%     sig_defect_pO2_sweep_cell = cell(num_rows, num_columns); % Pre-allocate correctly
%     % Populate the first row with the headers
%     sig_defect_pO2_sweep_cell(1,:) = cellstr(sig_defect_headers);
% end
% 
% % Populate the rest of the rows
% sig_defect_pO2_sweep_cell(pIdx + 1, :) = num2cell(sig_defects_out); % Assign data row-wise

%% Storing the Pressure Sweeps for all chargestates

% Creating a cell to store each iteration of pressure sweep
if pIdx == 1
    % Pre-allocate the cell
    num_rows = length(conditions.pO2_range) + 1; % Include the header row
    num_columns = length(all_chargestate_headers);   % Match number of headers
    all_chargestate_pO2_sweep_cell = cell(num_rows, num_columns); % Pre-allocate correctly
    % Populate the first row with the headers
    all_chargestate_pO2_sweep_cell(1,:) = cellstr(all_chargestate_headers);
end

% Populate the rest of the rows
all_chargestate_pO2_sweep_cell(pIdx + 1, :) = num2cell(all_chargestates_out); % Assign data row-wise

%% Storing the Pressure Sweeps for significant chargestates

% % Creating a cell to store each iteration of pressure sweep
% if pIdx == 1
%     % Pre-allocate the cell
%     num_rows = length(conditions.pO2_range) + 1; % Include the header row
%     num_columns = length(sig_chargestate_headers);   % Match number of headers
%     sig_chargestates_pO2_sweep_cell = cell(num_rows, num_columns); % Pre-allocate correctly
%     % Populate the first row with the headers
%     sig_chargestates_pO2_sweep_cell(1,:) = cellstr(sig_chargestate_headers);
% end
% 
% % Populate the rest of the rows
% sig_chargestates_pO2_sweep_cell(pIdx + 1, :) = num2cell(sig_chargestates_out); % Assign data row-wise
    

%% clean up memory
% clear equilib_save_fname
clear all_defects_out all_chargestates_out all_defects_headers all_chargestates_headers all_defects_cell all_chargstates_cell sig_defects_index sig_chargestates_index sig_defects_out sig_chargestates_out sig_defects_headers sig_chargestates_headers sig_defects_cell sig_chargestates_cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% end full equilibrium calc %%%

 end % end of pressure looping for loop










% % %% Task Trange 2: quenching from the annealing temperatures from example 1  %%%%%%%%%%%%%%%
% % 
% % %%% note this assumes equilib_dark_sol is still in memory - could also
% % %%% add a line here to load a prior solution from a mat file.
% % 
% % tic
% % [fullquench_dark_sol] = defect_fullquench_dark(equilib_dark_sol, conditions, defects);
% % toc
% % 
% % % these variables don't exist for quenching so set them to zero to keep the
% % % file format the same with the full equilibrium calc
% % fullquench_dark_sol.element_bal_err = zeros(size(fullquench_dark_sol.charge_bal_err));
% % fullquench_dark_sol.tot_bal_err = zeros(size(fullquench_dark_sol.charge_bal_err));
% % fullquench_dark_sol.mu = zeros(size(equilib_dark_sol.mu));
% % EvT_fullquench_vec = conditions.EvT_fullquench*ones(size(fullquench_dark_sol.charge_bal_err));
% % EcT_fullquench_vec = conditions.EcT_fullquench*ones(size(fullquench_dark_sol.charge_bal_err));
% % 
% % %%%%%%%%%%% save the outputs for quenching.  Yes some of these are the
% % %%%%%%%%%%% same as computed for the equilibrium solution but for sake
% % %%%%%%%%%%% of debugging/proofreading were jsut doing it here again
% % all_defects_out = cat(2,fullquench_dark_sol.T_equilibrium, fullquench_dark_sol.charge_bal_err, fullquench_dark_sol.element_bal_err, fullquench_dark_sol.tot_bal_err, fullquench_dark_sol.mu, conditions.EvT_equilibrium', conditions.EcT_equilibrium', fullquench_dark_sol.EFn, fullquench_dark_sol.EFp, fullquench_dark_sol.Nd, fullquench_dark_sol.Na, fullquench_dark_sol.n, fullquench_dark_sol.p, fullquench_dark_sol.sth1, fullquench_dark_sol.sth2, fullquench_dark_sol.defects);
% % all_chargestates_out = cat(2,fullquench_dark_sol.T_equilibrium, fullquench_dark_sol.charge_bal_err, fullquench_dark_sol.element_bal_err, fullquench_dark_sol.tot_bal_err, fullquench_dark_sol.mu, conditions.EvT_equilibrium', conditions.EcT_equilibrium', fullquench_dark_sol.EFn, fullquench_dark_sol.EFp, fullquench_dark_sol.Nd, fullquench_dark_sol.Na, fullquench_dark_sol.n, fullquench_dark_sol.p, fullquench_dark_sol.sth1, fullquench_dark_sol.sth2, fullquench_dark_sol.chargestates);
% % all_defect_headers = ['T_equilibrium (K)' 'charge_bal_err' 'element_bal_err' 'tot_bal_err' defects.mu_names_with_units 'Ev(T) (eV)' 'Ec(T) (eV)' 'Efn(eV)' 'Efp (eV)' 'Nd (#/cm3)' 'Na (#/cm3)' 'n (#/cm3)' 'p (#/cm3)' 'sth1 (#/cm3)' 'sth2 (#/cm3)' defects.defect_names'];
% % all_chargestate_headers = ['T_equilibrium (K)' 'charge_bal_err' 'element_bal_err' 'tot_bal_err' defects.mu_names_with_units 'Ev(T) (eV)' 'Ec(T) (eV)' 'Efn(eV)' 'Efp (eV)' 'Nd (#/cm3)' 'Na (#/cm3)' 'n (#/cm3)' 'p (#/cm3)' 'sth1 (#/cm3)' 'sth2 (#/cm3)' defects.chargestate_names'];
% % all_defects_cell = [cellstr(all_defect_headers); num2cell(all_defects_out)];
% % all_chargestates_cell = [cellstr(all_chargestate_headers); num2cell(all_chargestates_out)];
% % 
% % % throw out all the ones with insignificant concentrations and save
% % sig_defects_index = find(max(fullquench_dark_sol.defects,[],1)>=conditions.save_min);
% % sig_chargestates_index = find(max(fullquench_dark_sol.chargestates,[],1)>=conditions.save_min);
% % sig_defects_out = cat(2,fullquench_dark_sol.T_equilibrium, fullquench_dark_sol.charge_bal_err, fullquench_dark_sol.element_bal_err, fullquench_dark_sol.tot_bal_err, fullquench_dark_sol.mu, conditions.EvT_equilibrium', conditions.EcT_equilibrium', fullquench_dark_sol.EFn, fullquench_dark_sol.EFp, fullquench_dark_sol.Nd, fullquench_dark_sol.Na, fullquench_dark_sol.n, fullquench_dark_sol.p, fullquench_dark_sol.sth1, fullquench_dark_sol.sth2, fullquench_dark_sol.defects(:,sig_defects_index));
% % sig_chargestates_out = cat(2,fullquench_dark_sol.T_equilibrium, fullquench_dark_sol.charge_bal_err, fullquench_dark_sol.element_bal_err, fullquench_dark_sol.tot_bal_err, fullquench_dark_sol.mu, conditions.EvT_equilibrium', conditions.EcT_equilibrium', fullquench_dark_sol.EFn, fullquench_dark_sol.EFp, fullquench_dark_sol.Nd, fullquench_dark_sol.Na, fullquench_dark_sol.n, fullquench_dark_sol.p, fullquench_dark_sol.sth1, fullquench_dark_sol.sth2, fullquench_dark_sol.chargestates(:,sig_chargestates_index));
% % sig_defects_headers = ['T_equilibrium (K)' 'charge_bal_err' 'element_bal_err' 'tot_bal_err' defects.mu_names_with_units 'Ev(T) (eV)' 'Ec(T) (eV)' 'Efn(eV)' 'Efp (eV)' 'Nd (#/cm3)' 'Na (#/cm3)' 'n (#/cm3)' 'p (#/cm3)' 'sth1 (#/cm3)' 'sth2 (#/cm3)' defects.defect_names(sig_defects_index)'];
% % sig_chargestate_headers = ['T_equilibrium (K)' 'charge_bal_err' 'element_bal_err' 'tot_bal_err' defects.mu_names_with_units 'Ev(T) (eV)' 'Ec(T) (eV)' 'Efn(eV)' 'Efp (eV)' 'Nd (#/cm3)' 'Na (#/cm3)' 'n (#/cm3)' 'p (#/cm3)' 'sth1 (#/cm3)' 'sth2 (#/cm3)' defects.chargestate_names(sig_chargestates_index)'];
% % sig_defects_cell = [cellstr(sig_defect_headers); num2cell(sig_defects_out)];
% % sig_chargestates_cell = [cellstr(sig_chargestate_headers); num2cell(sig_chargestates_out)];
% % 
% % % write to a tab delim text file
% % % writematrix(equilib_dark_sol_out,strcat(conditions.save_pname,conditions.equilib_save_fname),'Delimiter','tab','WriteMode','overwrite')
% % if strcmp(conditions.save_files_flag,'All')
% %     writecell(all_defects_cell,strcat(conditions.save_pname,conditions.fullquench_save_fname,'_all_defects'),'Delimiter','tab','WriteMode','overwrite');  % write results with column headers
% %     writecell(all_chargestates_cell,strcat(conditions.save_pname,conditions.fullquench_save_fname,'_all_chargestates'),'Delimiter','tab','WriteMode','overwrite');
% % elseif strcmp(conditions.save_files_flag,'Sig_Only')
% %     writecell(sig_defects_cell,strcat(conditions.save_pname,conditions.fullquench_save_fname,'_sig_defects'),'Delimiter','tab','WriteMode','overwrite');
% %     writecell(sig_chargestates_cell,strcat(conditions.save_pname,conditions.fullquench_save_fname,'_sig_chargestates'),'Delimiter','tab','WriteMode','overwrite');
% % end
% % disp('Wrote quenched dark solution as text files')
% % 
% % % write the whole solution structure to a mat file
% % save(strcat(conditions.save_pname,conditions.fullquench_save_fname))
% % disp('Wrote quenched dark solution to .mat file')
% % 
% % % clean up memory
% % % clear fullquench_save_fname
% % clear  all_defects_out all_chargestates_out all_defects_headers all_chargestates_headers all_defects_cell all_chargstates_cell sig_defects_index sig_chargestates_index sig_defects_out sig_chargestates_out sig_defects_headers sig_chargestates_headers sig_defects_cell sig_chargestates_cell  EvT_fullquench_vec EcT_fullquench_vec
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%% end quenched calc %%%%%%%%%%%%%%%%%
% 
% 
% 
% 
% %% grouping defects into categories.  
% %% in the .mat file, each row is a list of column numbers from the main defect database to sum into a group 
% if strcmp(conditions.defect_group_flag,'On')
%     temp = load(strcat(conditions.defect_grouping_pname,conditions.defect_grouping_fname));  % this loads 3 variables into a structue "temp'
%     conditions.group_together_columns = temp.group_together_columns;
%     conditions.group_names = temp.group_names;
%     conditions.num_groups = temp.num_groups;
%     grouped_defects = zeros(numel(equilib_dark_sol.T_equilibrium), conditions.num_groups);
% 
%     for index = 1:conditions.num_groups
%         cols_to_sum = conditions.group_together_columns(index, conditions.group_together_columns(index, :) > 0);
%         grouped_defects(:, index) = sum(equilib_dark_sol.defects(:, cols_to_sum), 2);
%     end
% 
%     equilib_grouped_STH = equilib_dark_sol.sth1 + equilib_dark_sol.sth2;
%     quenched_grouped_STH = fullquench_dark_sol.sth1 + fullquench_dark_sol.sth2;
% 
%     grouping_headers = ['T_equilibrium (K)' 'n equilib (#/cm3)' 'p equilib (#/cm3)' 'STH equilib (#/cm3)' 'n quench (#/cm3)' 'p quench (#/cm3)' 'STH quench (#/cm3)' conditions.group_names'];
%     grouping_defects_out = cat(2,equilib_dark_sol.T_equilibrium, equilib_dark_sol.n, equilib_dark_sol.p, equilib_grouped_STH, fullquench_dark_sol.n, fullquench_dark_sol.p,  quenched_grouped_STH, grouped_defects);
%     grouping_defects_cell = [cellstr(grouping_headers); num2cell(grouping_defects_out)];
% 
%     writecell(grouping_defects_cell,strcat(conditions.save_pname, conditions.grouping_save_fname),'Delimiter','tab','WriteMode','overwrite');  % write results with column headers
% 
%     clear temp grouping_save_fname grouping_defects_cell grouping_defects_out grouping_headers grouped_defects index cols_to_sum grouped_STH
% end
% 
% 
% 
% 
% 
% %% Display results
% 
% 
% %% plot total defect concentrations - these will be same for equilib and quenched except for n and p
% figure(2)
% f2line(1) = plot(conditions.T_equilibrium,log10(equilib_dark_sol.n),'o','DisplayName','n equilib','MarkerEdgeColor',"r",'MarkerFaceColor',"r");  % full circles for equilib, open for fullquench
% f2line(2) = plot(conditions.T_equilibrium,log10(equilib_dark_sol.p),'o','DisplayName','p equilib','MarkerEdgeColor',"b",'MarkerFaceColor',"b");
% f2line(3) = plot(conditions.T_equilibrium,log10(equilib_dark_sol.sth1),'o','DisplayName','sth1 equilib','MarkerEdgeColor',"g",'MarkerFaceColor',"g");
% f2line(4) = plot(conditions.T_equilibrium,log10(equilib_dark_sol.sth2),'o','DisplayName','sth2 equilib','MarkerEdgeColor',"g",'MarkerFaceColor',"g");
% f2line(5) = plot(conditions.T_equilibrium,log10(abs(equilib_dark_sol.Nd - equilib_dark_sol.Na)),'k-','DisplayName','|Nd-Na|');
% f2line(6) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.n),'o','DisplayName','n quench','MarkerEdgeColor',"r");
% f2line(7) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.p),'o','DisplayName','p quench','MarkerEdgeColor',"b");
% f2line(8) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.sth1),'o','DisplayName','sth1 quench','MarkerEdgeColor',"g");
% f2line(9) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.sth2),'o','DisplayName','sth2 quench','MarkerEdgeColor',"g");
% 
% for i = 1:9
%     f2line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f2line(i).DisplayName},size(f2line(i).XData)));
% end
% legend('n equilib','p equilib', 'sth1 equilib', 'sth2 equilib', '|Nd-Na|', 'n quench', 'p quench',  'sth1 quench', 'sth2 quench' )
% 
% 
% for i=1:defects.num_defects
%     if max(log10(equilib_dark_sol.defects(:,i))) >= conditions.plot_log10_min   % only plot defects that at some point rise above 1e10 /cm3 - this will cut down the legend size
%         f2line(i+5) = plot(conditions.T_equilibrium,log10(equilib_dark_sol.defects(:,i)),'DisplayName', defects.defect_names(i));
%         f2line(i+5).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f2line(i+5).DisplayName},size(f2line(i+5).XData)));
%     end
% end
% legend show
% title('Concentrations at Equilibrium Temperature K')
% xlabel('Equilibrium T (K)')
% ylabel('Log_{10} Concentrations (#/cm3)')
% ylim([conditions.plot_log10_min conditions.plot_log10_max])
% xlim([min(conditions.T_equilibrium) max(conditions.T_equilibrium)])
% datacursormode.Enable='on';   %turn on the datatip to give name of defect on mouse hover
% 
% 
% 
% 
% %% plot all chargestates for equilibrium
% figure(3)
% f3line(1) = plot(conditions.T_equilibrium,log10(equilib_dark_sol.n),'o','DisplayName','n equilib','MarkerEdgeColor',"r",'MarkerFaceColor',"r");  % full circles for equilib, open for fullquench
% f3line(2) = plot(conditions.T_equilibrium,log10(equilib_dark_sol.p),'o','DisplayName','p equilib','MarkerEdgeColor',"b",'MarkerFaceColor',"b");
% f3line(3) = plot(conditions.T_equilibrium,log10(equilib_dark_sol.sth1),'o','DisplayName','sth1 equilib','MarkerEdgeColor',"g",'MarkerFaceColor',"g");
% f3line(4) = plot(conditions.T_equilibrium,log10(equilib_dark_sol.sth2),'o','DisplayName','sth2 equilib','MarkerEdgeColor',"g",'MarkerFaceColor',"g");
% f3line(5) = plot(conditions.T_equilibrium,log10(abs(equilib_dark_sol.Nd - equilib_dark_sol.Na)),'k-','DisplayName','|Nd-Na|');
% f3line(6) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.n),'o','DisplayName','n quench','MarkerEdgeColor',"r");
% f3line(7) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.p),'o','DisplayName','p quench','MarkerEdgeColor',"b");
% f3line(8) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.sth1),'o','DisplayName','sth1 quench','MarkerEdgeColor',"g");
% f3line(9) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.sth2),'o','DisplayName','sth2 quench','MarkerEdgeColor',"g");
% 
% for i = 1:9
%     f3line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f3line(i).DisplayName},size(f3line(i).XData)));
% end
% legend('n equilib','p equilib', 'sth1 equilib', 'sth2 equilib', '|Nd-Na|', 'n quench', 'p quench',  'sth1 quench', 'sth2 quench' )
% 
% 
% for i=1:defects.num_chargestates
%     if max(log10(equilib_dark_sol.chargestates(:,i))) >= conditions.plot_log10_min   % only plot defects that at some point rise above 1e10 /cm3 - this will cut down the legend size
%         f3line(i+3) = plot(conditions.T_equilibrium,log10(equilib_dark_sol.chargestates(:,i)),'DisplayName', defects.chargestate_names(i));
%         f3line(i+3).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f3line(i+3).DisplayName},size(f3line(i+3).XData)));
%     end
% end
% legend show
% title('Concentrations at Equilibrium Temperature K')
% xlabel('Equilibrium T (K)')
% ylabel('Log_{10} Concentrations (#/cm3)')
% ylim([conditions.plot_log10_min conditions.plot_log10_max])
% xlim([min(conditions.T_equilibrium) max(conditions.T_equilibrium)])
% datacursormode.Enable='on';   %turn on the datatip to give name of defect on mouse hover
% 
% 
% % 
% % % plot all chargestates for quenching - chargestates get changed for quenching
% % figure(4)
% % f4line(1) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.n),'o','DisplayName','n','MarkerEdgeColor',"r");
% % f4line(2) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.p),'o','DisplayName','p','MarkerEdgeColor',"b");
% % f4line(3) = plot(conditions.T_equilibrium,log10(abs(fullquench_dark_sol.Nd - fullquench_dark_sol.Na)),'k-','DisplayName','|Nd-Na|');
% % for i = 1:3
% %     f4line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f4line(i).DisplayName},size(f4line(i).XData)));
% % end
% % legend('n','p','|Nd-Na|')
% % for i=1:defects.num_chargestates
% %     if max(log10(fullquench_dark_sol.chargestates(:,i))) >= conditions.plot_log10_min   % only plot defects that at some point rise above 1e10 /cm3 - this will cut down the legend size
% %         f4line(i+3) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.chargestates(:,i)),'DisplayName', defects.chargestate_names(i));
% %         f4line(i+3).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f4line(i+3).DisplayName},size(f4line(i+3).XData)));
% %     end
% % end
% % legend show
% % title('Concentrations at Equilibrium Temperature K')
% % xlabel('Equilibrium T (K)')
% % ylabel('Log_{10} Concentrations (#/cm3)')
% % ylim([conditions.plot_log10_min conditions.plot_log10_max])
% % xlim([min(conditions.T_equilibrium) max(conditions.T_equilibrium)])
% % datacursormode.Enable='on';   %turn on the datatip to give name of defect on mouse hover
% 
% 
% 
% %% Plot the  band diagram vs T Ec, Ev, EFn and EFp vs T
% figure(5)
% 
% plot(conditions.T_equilibrium,equilib_dark_sol.EFp,'b-')   % note in dark Efn and Efp will lie over eachother
% plot(conditions.T_equilibrium,equilib_dark_sol.EFn,'r-')
% % plot(conditions.T_equilibrium,fullquench_dark_sol.EFp,'--b')
% % plot(conditions.T_equilibrium,fullquench_dark_sol.EFn,'--r')
% 
% %% handle the T dependent band edges or constant ones
% % if strcmp(conditions.T_dep_bands_flag,'On')
% plot(conditions.T_equilibrium,conditions.EcT_equilibrium,'k-')
% plot(conditions.T_equilibrium,conditions.EvT_equilibrium,'k-')
% % elseif strcmp(conditions.T_dep_bands_flag,'Off')
% plot(conditions.T_equilibrium,conditions.EcRef*ones(size(conditions.T_equilibrium)),'k--')
% plot(conditions.T_equilibrium,conditions.EvRef*ones(size(conditions.T_equilibrium)),'k--')
% % else
% % error('conditions.T_dep_bands_flag must be on or off')
% % end
% 
% title('Fermi Levels for Equilibrium and Quenching vs T')
% legend('EFp(T_{Equilib})','EFn(T_{Equilib})', 'EFp(T_{Quench})', 'EFn(T_{Quench})','Ec(T_{Equilib})','Ev(T_{Equilib})','Ec(T_{Quench})','Ev(T_{Quench})')
% xlabel('Equilibrium T (K)')
% ylabel('Energies (eV)')
% 
% 
% 
% 
% 
% %% Plot the mu values for the solution
% figure(6)
% 
% %% handle the T dependent band edges or constant ones
% % for i=1:defects.numelements
% %     plot(conditions.T_equilibrium, equilib_dark_sol.mu)
% % end
% 
% f6line(1) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,1),'DisplayName','Ga');
% f6line(2) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,2),'DisplayName','O');
% f6line(3) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,3),'DisplayName','Si');
% f6line(4) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,4),'DisplayName','H');
% f6line(5) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,5),'DisplayName','Fe');
% f6line(6) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,6),'DisplayName','Sn');
% f6line(7) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,7),'DisplayName','Cr');
% f6line(8) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,8),'DisplayName','Ti');
% f6line(9) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,9),'DisplayName','Ir');
% f6line(10) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,10),'DisplayName','Mg');
% f6line(11) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,11),'DisplayName','Ca');
% f6line(12) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,12),'DisplayName','Zn');
% f6line(13) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,13),'DisplayName','Co');
% f6line(14) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,14),'DisplayName','Zr');
% f6line(15) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,15),'DisplayName','Hf');
% f6line(16) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,16),'DisplayName','Ta');
% f6line(17) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,17),'DisplayName','Ge');
% 
% for i = 1:17
%     f6line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f6line(i).DisplayName},size(f6line(i).XData)));
% end
% title('Chemical Potentials for Elements at Equilibrium and Quenching vs T')
% legend('\mu_{Ga}', '\mu_{O}', '\mu_{Si}', '\mu_{H}', '\mu_{Fe}', '\mu_{Sn}', '\mu_{Cr}', '\mu_{Ti}', '\mu_{Ir}', '\mu_{Mg}', '\mu_{Ca}', '\mu_{Zn}', '\mu_{Co}', '\mu_{Zr}', '\mu_{Hf}', '\mu_{Ta}', '\mu_{Ge}');
% xlabel('Equilibrium T (K)')
% ylabel('Chem Potentials (eV)')
% ylim([-12 0])
% 
% 
% 
% 
% %% plot the stoichiometry from equilibrium
% figure(7)
% 
% for i = 1:conditions.num_elements
%     f7line(i) = plot(conditions.T_equilibrium,log10(equilib_dark_sol.element_totals(:,i)),'DisplayName',defects.elementnames(i));
%     f7line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f7line(i).DisplayName},size(f7line(i).XData)));
% end
% 
% %
% % f7line(2) = plot(conditions.T_equilibrium,log10(num_O),'DisplayName','O');
% % f7line(3) = plot(conditions.T_equilibrium,log10(num_Si),'DisplayName','Si');
% % f7line(4) = plot(conditions.T_equilibrium,log10(num_H),'DisplayName','H');
% % f7line(5) = plot(conditions.T_equilibrium,log10(num_Fe),'DisplayName','Fe');
% % f7line(6) = plot(conditions.T_equilibrium,log10(num_Sn),'DisplayName','Sn');
% % f7line(7) = plot(conditions.T_equilibrium,log10(num_Cr),'DisplayName','Cr');
% % f7line(8) = plot(conditions.T_equilibrium,log10(num_Ti),'DisplayName','Ti');
% % f7line(9) = plot(conditions.T_equilibrium,log10(num_Ir),'DisplayName','Ir');
% % f7line(10) = plot(conditions.T_equilibrium,log10(num_Mg),'DisplayName','Mg');
% % f7line(11) = plot(conditions.T_equilibrium,log10(num_Ca),'DisplayName','Ca');
% % f7line(12) = plot(conditions.T_equilibrium,log10(num_Zn),'DisplayName','Zn');
% % f7line(13) = plot(conditions.T_equilibrium,log10(num_Co),'DisplayName','Co');
% % f7line(14) = plot(conditions.T_equilibrium,log10(num_Zr),'DisplayName','Zr');
% % f7line(15) = plot(conditions.T_equilibrium,log10(num_Hf),'DisplayName','Hf');
% % f7line(16) = plot(conditions.T_equilibrium,log10(num_Ta),'DisplayName','Ta');
% % f7line(17) = plot(conditions.T_equilibrium,log10(num_Ge),'DisplayName','Ge');
% 
% 
% title('Atomic Concentrations')
% legend('Ga','O','Si','H','Fe','Sn','Cr', 'Ti', 'Ir', 'Mg', 'Ca', 'Zn', 'Co', 'Zr', 'Hf', 'Ta', 'Ge')
% xlabel('Equilibrium T (K)')
% ylabel('Log_{10} Concentrations (#/cm3)')
% ylim([conditions.plot_log10_min conditions.plot_log10_max])
% 
% 
% 
% 
% 
% % analyze the distribution of each element across defects
% 
% % here just using it for the frozen ones but could extend this for all
% % the elements inclduing those with fixed mu
% 
% for ii = 1:numel(conditions.indices_of_fixed_elements)
% 
%     figure(7+ii)
%     clf
%     hold on
%     [defects_with_element_fraction,~, defects_with_element_names, ~] = contains_element(equilib_dark_sol, defects, conditions.indices_of_fixed_elements(ii));
% 
%     for i=1:size(defects_with_element_fraction,2)
%         plot(equilib_dark_sol.T_equilibrium,defects_with_element_fraction(:,i))
%         %     eval(strcat("f",int2str(7+ii),"line(i) = plot(equilib_dark_sol.T_equilibrium,defects_with_element_fraction(:,i),'DisplayName',defects_with_element_names(i));"))
%         %     eval(strcat("f",int2str(7+ii),"line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f",int2str(7+ii),"line(i).DisplayName},size(f"int2str(7+ii),"line(i).XData)));"));
%     end
% 
%     title('Fraction of Element in Defects')
%     legend(defects_with_element_names')
%     xlabel('Equilibrium T (K)')
%     ylabel('Fraction of Element in Each Defect')
%     ylim([0 1])
% 
% end
% 
% 
% 
% 
% 
% %%% clean up memory
% clear datacursormode f1line f11line f2line f3line f4line f5 line f6line f7line i
% clear num_Ga num_O num_Si num_H num_Fe num_Sn num_Cr num_Ti num_Ir num_Mg num_Ca num_Zn num_Co num_Zr num_Hf num_Ta num_Ge ii defects_with_element_fraction defects_with_element_names
% clear all_elements_out all_defects_cell all_chargestates_out all_defects_out all_defects_cell all_defect_headers all_elements_headers f4line all_elements_cell sig_chargestates_cell sig_chargestates_out element_totals_cell element_totals_headers elements_save_fname
% clear all_chargestate_headers all_chargestates_cell save_pname sig_chargestate_headers sig_defect_headers
% %%%%%%%%%%%%% end Tasks 1 and 2 (equilib and fullquenching)
% 
% 
% 
% % end % end for loop over Si doping
% 
% % end  % end the for loop with index "index"
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % % %%%%%%%%%%%%%%% Task Trange 3 - quenching but setting interstitials generated at the equilibrium temp to zero   %%%%%%%%%%%%%%%
% % %
% % % % take out the interstitials by hand.  These defect and chargestate numbers
% % % % only work for the full Varley native defect set (not validation)
% % %
% % % equilibrium_minus_interstitials_dark_sol = equilib_dark_sol;
% % %
% % % if strcmp(defect_set_flag, 'validation')
% % %     %% If using the validation data set the interstitials are column 2 and 4.
% % %     equilibrium_minus_interstitials_dark_sol.defects(:,2)=0;    % Ga_i
% % %     equilibrium_minus_interstitials_dark_sol.chargestates(:,6:8)=0;
% % %     equilibrium_minus_interstitials_dark_sol.defects(:,4)=0;   % O_i
% % %     equilibrium_minus_interstitials_dark_sol.chargestates(:,12:14)=0;
% % % elseif strcmp(defect_set_flag, 'full_Varley_set')
% % %     %% if using the full Varley set, the interstitials are defects numbers 6, 13, and 14, and charge states 26,27,28 (Gai), and 52-56 (Oi and Osi)
% % %     equilibrium_minus_interstitials_dark_sol.defects(:,6)=0;    % Ga_i
% % %     equilibrium_minus_interstitials_dark_sol.chargestates(:,26:28)=0;
% % %     equilibrium_minus_interstitials_dark_sol.defects(:,13)=0;   % O_i
% % %     equilibrium_minus_interstitials_dark_sol.chargestates(:,52:54)=0;
% % %     equilibrium_minus_interstitials_dark_sol.defects(:,14)=0;   % O_si
% % %     equilibrium_minus_interstitials_dark_sol.chargestates(:,55:56)=0;
% % % end
% % %
% % % tic
% % % [quench_minus_interstitials_dark_sol] = defect_fullquench_dark(equilibrium_minus_interstitials_dark_sol, conditions, defects);
% % % toc
% % %
% % % quench_minus_interstitials_dark_sol_out = cat(2,quench_minus_interstitials_dark_sol.T_equilibrium, quench_minus_interstitials_dark_sol.Nd, quench_minus_interstitials_dark_sol.Na, quench_minus_interstitials_dark_sol.EFn, quench_minus_interstitials_dark_sol.EFp, quench_minus_interstitials_dark_sol.n, quench_minus_interstitials_dark_sol.p, quench_minus_interstitials_dark_sol.defects);
% % % save_name = strcat(conditions.name_root,'_quench_minus_interstitials_sol.dat');
% % % writematrix(quench_minus_interstitials_dark_sol_out,save_name,'Delimiter','tab','WriteMode','overwrite')
% % % disp(strcat(save_name,' written.'))
% % %
% % %
% % % %% plot total defect concentrations
% % % figure(6)
% % % f6line(1) = plot(conditions.T_equilibrium,log10(quench_minus_interstitials_dark_sol.n),'ro');
% % % f6line(2) = plot(conditions.T_equilibrium,log10(quench_minus_interstitials_dark_sol.p),'bo');
% % % f6line(3) = plot(conditions.T_equilibrium,log10(abs(quench_minus_interstitials_dark_sol.Nd - quench_minus_interstitials_dark_sol.Na)),'k-');
% % % for i = 1:3
% % %     f6line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f6line(i).DisplayName},size(f6line(i).XData)));
% % % end
% % % legend('n','p','|Nd-Na|')
% % % for i=1:defects.num_defects
% % %     f6line(i+3) = plot(conditions.T_equilibrium,log10(quench_minus_interstitials_dark_sol.defects(:,i)),'DisplayName', defects.defect_names(i));
% % %     f6line(i+3).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f6line(i+3).DisplayName},size(f6line(i+3).XData)));
% % % end
% % % legend show
% % % title('Concentrations for Quenching all but O_i and Ga_i to 300 K')
% % % xlabel('T for Equilibriation (K)')
% % % ylabel('Log_{10} Concentrations at 300 K (#/cm3)')
% % % ylim([10 23])
% % % datacursormode.Enable='on';   %turn on the datatip to give name of defect on mouse hover
% % %
% % %
% % %
% % % %% add quenching minus interstitials results to Fig. 4
% % % figure(4)
% % % plot(conditions.T_equilibrium,quench_minus_interstitials_dark_sol.EFp,'k-')
% % % plot(conditions.T_equilibrium,quench_minus_interstitials_dark_sol.EFn,'g-')
% % % legend('Ec(T)','Ec(T_{quench})','Ev(T)','EFp(T) Equilib','EFn(T) Equilib', 'EFp(T) Quench', 'EFn(T) Quench', 'EFp(T) Quench, no interstitials', 'EFn(T) Quench, no interstitials')
% % % xlabel('Equilibrium T (K)')
% % % ylabel('Ec, Ev, EFn, EFp (eV)')
% %
% %
% %
% %
% % %%%%% the tasks below here reqire a single temperature
% % if size(conditions.T_equilibrium,2)~=1
% %     disp('calculation of defects vs doping at fixed T requires only one T be in conditions.T_equilibrium - chose one T rather than a vector of T values')
% %     yn = input('do you want to run the other tasks at a fixed temperature (must enter string y or n)');
% %
% %     if isempty(yn)
% %         error('stopping... other parameter sweep tasks will be skipped')
% %
% %     elseif strcmp(yn,'y')
% %         conditions.T_equilibrium = input('Enter one fixed temperature (K) for the rest of the calculations.  default T is set to 1300 K if you hit enter');
% %         if isempty(conditions.T_equilibrium)
% %             conditions.T_equilibrium = 1300;
% %             for i=1:10
% %                 disp('T set to 1000 C or 1300 K automatically')
% %             end
% %         else
% %             disp('calculation of defects vs doping at fixed T requires only one T be in conditions.T_equilibrium - chose one T rather than a vector of T values')
% %         end
% %     else
% %         error('stopping... other parameter sweep tasks will be skipped')
% %     end
% % end
% %
% %
% %
% %
% % %%% Task Doping range 1:  explore defects vs doping at constant T %%%%%
% %
% %
% % % make a temp holder for the shallow doping set at the start of the script,
% % % reset back to that value once done.
% % temp_Nd = conditions.Nd;
% % temp_Na = conditions.Na;
% %
% % log_doping = 10:0.1:20;  %log10 of the doping range you want to sweep over
% % num_dope = size(log_doping,2);
% % defect_holder = zeros(num_dope,defects.num_defects);
% % % chargestate_holder = zeros(num_dope,defects.num_chargestates);
% % atom_holder = zeros(num_dope,6);    % change the 5 to a 6 once Sn added
% % equilib_carrier_holder = zeros(num_dope,3);
% % quench_carrier_holder = zeros(num_dope,3);
% %
% % for i = 1:num_dope
% %     conditions.Nd = 10^log_doping(i);  %vary n-type shallow doping
% % %     conditions.Na = 10^log_doping(i);  %vary p-type shallow doping
% %     [equilib_dark_sol] = defect_equilibrium_dark(conditions, defects);
% %     defect_holder(i,:) = equilib_dark_sol.defects;
% % %     chargestate_holder(i,:) = equilib_dark_sol.defects;
% %     equilib_carrier_holder(i,1) = equilib_dark_sol.n;
% %     equilib_carrier_holder(i,2) = equilib_dark_sol.p;
% %     equilib_carrier_holder(i,3) = abs(equilib_dark_sol.Nd - equilib_dark_sol.Na);  %no need to duplicate this for quenching
% %     [fullquench_dark_sol] = defect_fullquench_dark(equilib_dark_sol, conditions, defects);
% %     quench_carrier_holder(i,1) = fullquench_dark_sol.n;
% %     quench_carrier_holder(i,2) = fullquench_dark_sol.p;
% %     [Ga_O_stoich, atom_holder(i,1), atom_holder(i,2), atom_holder(i,3), atom_holder(i,4), atom_holder(i,5), atom_holder(i,6)] =  Ga2O3_stoich(equilib_dark_sol, conditions, defects);
% % end
% %
% %
% % figure(7)
% % f7line(1) = plot(log_doping,log10(equilib_carrier_holder(:,1)),'ro','DisplayName','n_equilib');
% % f7line(2) = plot(log_doping,log10(equilib_carrier_holder(:,2)),'bo','DisplayName','p_equilib');
% % f7line(3) = plot(log_doping,log10(equilib_carrier_holder(:,3)),'k-','DisplayName','|Nd-Na|');  %this doesnt change for quenching so dont need to replot
% % f7line(4) = plot(log_doping,log10(quench_carrier_holder(:,1)),'rx','DisplayName','n_quech');
% % f7line(5) = plot(log_doping,log10(quench_carrier_holder(:,2)),'bx','DisplayName','p_quench');
% % % for i = 1:5
% % %     f7line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f7line(i).DisplayName},size(f7line(i).XData)));
% % % end
% % legend('n_{equilib}','p_{equilib}','|Nd-Na|','n_{quench}','p_{quench}')
% %
% % for i=1:defects.num_defects  %defects dont change for quenching so only need to plot once
% %     f7line(i+5) = plot(log_doping,log10(defect_holder(:,i)),'DisplayName', defects.defect_names(i));
% %     f7line(i+5).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f7line(i+5).DisplayName},size(f7line(i+5).XData)));   %the i+5 here is becasue we have 5 lines above for carrier densities
% % end
% % legend show
% % title('Concentrations at Fixed T vs Shallow Doping')
% % xlabel('Log_{10}(Nd) (#/cm^3)')
% % ylabel('Log_{10} Concentrations (#/cm^3)')
% % ylim([10 23])
% % datacursormode.Enable='on';   %turn on the datatip to give name of defect on mouse hover
% %
% %
% %
% % figure(8)
% % f8line(1) = plot(log_doping,log10(equilib_carrier_holder(:,1)),'ro','DisplayName','n_equilib');
% % f8line(2) = plot(log_doping,log10(equilib_carrier_holder(:,2)),'bo','DisplayName','p_equilib');
% % f8line(3) = plot(log_doping,log10(equilib_carrier_holder(:,3)),'k-','DisplayName','|Nd-Na|');  % this doesnt change for quenching so dont need to replot
% % f8line(4) = plot(log_doping,log10(quench_carrier_holder(:,1)),'rx','DisplayName','n_quech');
% % f8line(5) = plot(log_doping,log10(quench_carrier_holder(:,2)),'bx','DisplayName','p_quench');
% % for i = 1:5
% %     f8line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f8line(i).DisplayName},size(f8line(i).XData)));
% % end
% % legend('n_{equilib}','p_{equilib}','|Nd-Na|','n_{quench}','p_{quench}')
% %
% % for i=1:6  % defects dont change for quenching so only need to plot once.  Plot the atomic concentrations
% %     f8line(i+5) = plot(log_doping,log10(atom_holder(:,i)),'DisplayName', defects.elementnames(i));
% %     f8line(i+5).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f8line(i+5).DisplayName},size(f8line(i+5).XData)));   %the i+5 here is becasue we have 5 lines above for carrier densities
% % end
% % legend show
% % title('Concentrations at Fixed T vs Shallow Doping')
% % xlabel('log_{10}(Nd) (#/cm^3)')
% % ylabel('Log_{10} Concentrations (#/cm^3)')
% % ylim([10 23])
% % datacursormode.Enable='on';   %turn on the datatip to give name of defect on mouse hover
% %
% % %% reset the Nd and Na now that we're done
% % conditions.Nd = temp_Nd;
% % conditions.Na = temp_Na;
% % clear temp_Na temp_Nd f8line f7line datacursormode
% % %%%%%%%%%%%%%%% end of Task 4 loop over doping
% % %%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %
% %
% %
% % % % %  Task 5 loop over mu for some element
% % %% to use this, you must select the right line in the for loop below for the selected element
% %
% % for i=1:10
% %     disp('dont forget to reset the conditions.muT_equilibrium values after running task 5')
% % end
% %
% % %%%%%%%%% take the mu set at top of script, reset it back after this task
% % temp_mu = conditions.muT_equilibrium;
% %
% % if size(conditions.T_equilibrium,2)~=1
% %     error('Task 5 calculation of defects vs mu at fixed T requires only one T be in conditions.T_equilibrium - chose one T rather than a vector of T values')
% % end
% %
% % mu_vec = -6:0.1:0;  % enter the range to search over
% % num_mu_in_mu_vec = size(mu_vec,2);
% % % chargestate_holder = zeros(num_mu_in_mu_vec,defects.num_chargestates);
% % defect_holder = zeros(num_mu_in_mu_vec,defects.num_defects);
% % atom_holder = zeros(num_mu_in_mu_vec,defects.numelements);    % 6 corresponds to the list Ga O Si H Fe Sn
% % equilib_carrier_holder = zeros(num_mu_in_mu_vec,3);
% % quench_carrier_holder = zeros(num_mu_in_mu_vec,3);
% %
% % for i = 1:num_mu_in_mu_vec
% %     conditions.muT_equilibrium(:,conditions.Task5_element_num_to_vary_mu) = mu_vec(i);   % assign the desired value of mu to the right element
% %     [equilib_dark_sol] = defect_equilibrium_dark(conditions, defects);  % run the calc
% %     defect_holder(i,:) = equilib_dark_sol.defects;
% % %     chargestate_holder(i,:) = equilib_dark_sol.defects;
% %     equilib_carrier_holder(i,1) = equilib_dark_sol.n;
% %     equilib_carrier_holder(i,2) = equilib_dark_sol.p;
% %     equilib_carrier_holder(i,3) = abs(equilib_dark_sol.Nd - equilib_dark_sol.Na);
% %     [fullquench_dark_sol] = defect_fullquench_dark(equilib_dark_sol, conditions, defects);
% %     quench_carrier_holder(i,1) = fullquench_dark_sol.n;
% %     quench_carrier_holder(i,2) = fullquench_dark_sol.p;
% %     quench_carrier_holder(i,3) = abs(fullquench_dark_sol.Nd - fullquench_dark_sol.Na);
% %     [Ga_O_stoich, atom_holder(i,1), atom_holder(i,2), atom_holder(i,3), atom_holder(i,4), atom_holder(i,5), atom_holder(i,6)] =  Ga2O3_stoich(equilib_dark_sol, conditions, defects);
% % end
% %
% % figure(9)
% % f9line(1) = plot(mu_vec,log10(equilib_carrier_holder(:,1)),'ro','DisplayName','n');
% % f9line(2) = plot(mu_vec,log10(equilib_carrier_holder(:,2)),'bo','DisplayName','p');
% % f9line(3) = plot(mu_vec,log10(equilib_carrier_holder(:,3)),'k-','DisplayName','|Nd-Na|');
% % f9line(4) = plot(mu_vec,log10(quench_carrier_holder(:,1)),'rx','DisplayName','n_quech');
% % f9line(5) = plot(mu_vec,log10(quench_carrier_holder(:,2)),'bx','DisplayName','p_quench');
% % for i = 1:5
% %     f9line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f9line(i).DisplayName},size(f9line(i).XData)));
% % end
% % legend('n_{equilib}','p_{equilib}','|Nd-Na|','n_{quench}','p_{quench}')
% % for i=1:defects.num_defects
% %     f9line(i+5) = plot(mu_vec,log10(defect_holder(:,i)),'DisplayName', defects.defect_names(i));
% %     f9line(i+5).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f9line(i+3).DisplayName},size(f9line(i+3).XData)));
% % end
% % legend show
% % title(strcat('Concentrations at Fixed T vs_',defects.elementnames(conditions.Task5_element_num_to_vary_mu),' Chemical Potential'))
% % xlabel(strcat(defects.elementnames(conditions.Task5_element_num_to_vary_mu),' Chemical potential (eV)'))
% % ylabel('Log_{10} Concentrations (#/cm^3)')
% % ylim([10 23])
% % datacursormode.Enable='on';   %turn on the datatip to give name of defect on mouse hover
% %
% %
% %
% %
% % figure(10)
% % f10line(1) = plot(mu_vec,log10(equilib_carrier_holder(:,1)),'ro','DisplayName','n_equilib');
% % f10line(2) = plot(mu_vec,log10(equilib_carrier_holder(:,2)),'bo','DisplayName','p_equilib');
% % f10line(3) = plot(mu_vec,log10(equilib_carrier_holder(:,3)),'k-','DisplayName','|Nd-Na|');  % this doesnt change for quenching so dont need to replot
% % f10line(4) = plot(mu_vec,log10(quench_carrier_holder(:,1)),'rx','DisplayName','n_quech');
% % f10line(5) = plot(mu_vec,log10(quench_carrier_holder(:,2)),'bx','DisplayName','p_quench');
% % for i = 1:5
% %     f10line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f10line(i).DisplayName},size(f10line(i).XData)));
% % end
% % legend('n_{equilib}','p_{equilib}','|Nd-Na|','n_{quench}','p_{quench}')
% %
% % for i=1:6  % defects dont change for quenching so only need to plot once.  Plot the atomic concentrations
% %     f10line(i+5) = plot(mu_vec,log10(atom_holder(:,i)),'DisplayName', defects.elementnames(i));
% %     f10line(i+5).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f10line(i+5).DisplayName},size(f10line(i+5).XData)));   %the i+5 here is becasue we have 5 lines above for carrier densities
% % end
% % legend show
% % title('Element Concentrations at Fixed T vs Shallow Doping')
% % xlabel(strcat(defects.elementnames(conditions.Task5_element_num_to_vary_mu),' Chemical potential (eV)'))
% % ylabel('Log_{10} Concentrations (#/cm^3)')
% % ylim([10 23])
% % datacursormode.Enable='on';   %turn on the datatip to give name of defect on mouse hover
% %
% % conditions.muT_equilibrium = temp_mu ;   % reset the mu values back to what they were before running task 5
% % clear temp_mu f9line f10line datacursormode
% %
% % %%%%%%%%%%%%%%%%%% end of Task 5 %%%%%%%%%%%%%%
% 
% 

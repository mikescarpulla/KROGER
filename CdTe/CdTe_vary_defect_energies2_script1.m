





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




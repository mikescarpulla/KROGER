%% clear workspace and figures   %%%%%%%%%%%%
clear all
clc


%% set Physical Constants  %%%%%%%%%
conditions.q = 1.602176565e-19;
conditions.h = 6.62606957e-34;
conditions.kB = 8.6173324e-5;
conditions.mo = 9.1093837e-31;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    conditions.defect_db_fname = 'CdTe_defects_Intuon_10062025.mat';
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




%%% here is a place you can on the fly modify any values in the saved
%%% database - for example to check the sensitivity of results on formation
%%% energies
%
% %%%%%%%%%%%%%%%%%% Manually change dHo for Vga and complexes - numbers based on Varley 031424 database
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

conditions.plot_log10_min = 12;   %%% set a threshold - if a quantity doesn't go above this at some temperature then dont bother plotting it
conditions.plot_log10_max = 20;   %%% set upper bound for plots too.  The real numbers are 10^x for these x's

% save to file options
% min cutoff concentration for significant concentrations
conditions.save_min = 1e12;  % set a minimum threshold for a defect or chargestate to be saved to the reduced output files (real number not log of it)

conditions.stoich_flag = 'CdTe';
% conditions.stoich_flag = 'Ga2O3';
% conditions.stoich_flag = 'GaN';

% save all the defects and chargestates or just the ones
% conditions.save_files_flag = 'All';
conditions.save_files_flag = 'Sig_Only';  % save only the ones above the save_min

% conditions.defect_group_flag = 'On';
conditions.defect_group_flag = 'Off';
if strcmp(conditions.defect_group_flag,'On')
    conditions.defect_grouping_fname = 'CdTe_grouping_defects.mat';
    % conditions.defect_grouping_fname = 'Ga2O3_grouping_defects.mat';
    conditions.defect_grouping_pname = strcat(pwd,'\');
end

% turn on Fermi-Dirac stats for bands or use Boltzmann stats
% conditions.Boltz_or_FD_flag = 'Boltz';
conditions.Boltz_or_FD_flag = 'FD';

% turn on or off site blocking statistics
% conditions.site_blocking_flag= 'On_Infinite';  % site blocking for infinite crystal
% conditions.site_blocking_flag= 'On_Finite';  % site blocking for finite crystal
conditions.site_blocking_flag= 'Off';

% turn on or off vibrational entropy.  Add 3kB per added atom, remove that for removed atoms.  G = H-TS so interstiitals should decrease G while vacacnies add it
% conditions.vib_ent_flag = '3kB';
conditions.vib_ent_flag = 'Quantum';
% conditions.vib_ent_flag = 'Off';

% Turn on or off T dependent band parameters
conditions.T_dep_bands_flag = 'On';
% conditions.T_dep_bands_flag = 'Off';

% give the relaxation energy for the STH
% conditions.sth_flag = 0; % this turns off the STH's
conditions.sth_flag = 0;  % this turns on the STH's
conditions.E_relax_sth1 = 0.0;  % relaxation energy for STH1's
conditions.E_relax_sth2 = 0.0;  % relaxation energy for STH2's

% set how to handle matrix chemical potentials
conditions.T_dep_matrix_mu_flag = 'On';      % using actual thermochemistry
% conditions.T_dep_matrix_mu_flag = 'Off';   % traditional way from most DFT papers - constant value vs T

% set how to handle defects with fixed concentrtations - either a
% constant value for all T or different values vs T
% conditions.T_dep_fixed_defect_flag = 'On';
conditions.T_dep_fixed_defect_flag = 'Off';

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
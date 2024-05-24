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


%% clear workspace and figures   %%%%%%%%%%%%
clear all
clc

conditions.figures_flag = 'On';
% conditions.figures_flag = 'Off';

if strcmp(conditions.figures_flag,'On')
    for i=1:7
        figure(i)
        clf
        hold on
    end
    clear i
end
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


[conditions.defect_db_fname,conditions.defect_db_pname] = uigetfile('*.mat','Select the defect database file you want to use');

if ~ischar(conditions.defect_db_fname)   %if the user doesnt pick one in the UIgetfile, do it manually here

    % defect_set_flag = 'Varley_validation';
    % defect_set_flag = 'Varley_validation_Sn';
    % defect_set_flag = 'Varley_native';
    % defect_set_flag = 'Varley_native_H';
    % defect_set_flag = 'Varley_native_Sn';
    % defect_set_flag = 'Varley_native_H_Sn';
    % defect_set_flag = 'Varley_all_defects';
    defect_set_flag = 'Varley_all_defects_new_032124';

    % load in the right set of defects from mat file
    if strcmp(defect_set_flag, 'Varley_validation')
        conditions.defect_set_fname = 'Ga2O3_Varley_validation_data.mat';
    elseif strcmp(defect_set_flag, 'Varley_validation_Sn')
        conditions.defect_set_fname = 'Ga2O3_Varley_validation_Sn_data.mat';
    elseif strcmp(defect_set_flag, 'Varley_native')
        conditions.defect_set_fname = 'Ga2O3_Varley_native_data.mat';
    elseif  strcmp(defect_set_flag, 'Varley_native_H')
        conditions.defect_set_fname = 'Ga2O3_Varley_native_H_data.mat';
    elseif  strcmp(defect_set_flag, 'Varley_native_Sn')
        conditions.defect_set_fname = 'Ga2O3_Varley_native_Sn_data_new.mat';
    elseif  strcmp(defect_set_flag, 'Varley_native_H_Sn')
        conditions.defect_set_fname = 'Ga2O3_Varley_native_H_Sn_data_new.mat';
    elseif  strcmp(defect_set_flag, 'Varley_all_defects')
        conditions.defect_set_fname = 'Ga2O3_Varley_all_defects.mat';
    elseif  strcmp(defect_set_flag, 'Varley_all_defects_new_032124')
        conditions.defect_set_fname = 'Ga2O3_Varley_all_defecs_new_032124.mat';
    else
        error('you have to use one of the databases with the right name')
    end
end

if ~ischar(conditions.defect_db_pname)
    conditions.defect_db_pname = strcat(pwd,'\');  % the pwd command leaves off the last \ so have to add it.
end

load(strcat(conditions.defect_db_pname,conditions.defect_db_fname))
defect_db_name_root = conditions.defect_db_fname(1:(size(conditions.defect_db_fname,2) - 4)    );  % this takes the filename and strips the . and extension - keep it until savenames established

if numel(fieldnames(defects)) < 16
    error('Defects database file seems to be missing some required parts - carefully check it and try again')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% choose folder and base filename to save calculation run outputs
[save_fname_root,save_pname] = uiputfile('*.dat','Select or type the filename and location for saving results');
save_fname_root_length = size(save_fname_root,2);

if ~ischar(save_fname_root)  %if user exits the UI window with no input
    save_fname_root = strcat(defect_db_name_root,strcat('_',strrep(char(datetime),' ','_')));  % default is to use the defect database name but append date and time to it
    save_pname = conditions.defect_db_pname;  %default is to use the same directory with the defect database in it
elseif strcmp(save_fname_root(save_fname_root_length-3),'.')
    save_fname_root = save_fname_root(1:(save_fname_root_length-4) );   % this takes the filename and strips the extension.  If you type a name without extension, the uiputfile adds .dat to it
    save_pname = conditions.defect_db_pname;  %default is to use the same directory with the defect database in it
end

conditions_save_fname = strcat(save_fname_root,'_calculation_conditions');
equilib_save_fname = strcat(save_fname_root,'_equilib_dark_sol');
fullquench_save_fname = strcat(save_fname_root,'_fullquench_dark_sol');
clear save_fname_root save_fname_root_length defect_set_flag defect_db_name_root
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% Options  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

conditions.plot_log10_min = 14;   %%% set a threshold - if a quantity doesn't go above this at some temperature then dont bother plotting it
conditions.plot_log10_max = 23;   %%% set upper bound for plots too.  The real numbers are 10^x for these x's

conditions.save_min = 1e14;  % set a minimum threshold for a defect or chargestate to be saved to the reduced output files

% turn on Fermi-Dirac stats for bands or use Boltzmann stats
% conditions.Boltz_or_FD_flag = 'Boltz';
conditions.Boltz_or_FD_flag = 'FD';

% turn on or off site blocking statistics
conditions.site_blocking_flag= 'On_Infinite';  % site blocking for infinite crystal
% conditions.site_blocking_flag= 'On_Finite';  % site blocking for finite crystal
% conditions.site_blocking_flag= 'Off';

% turn on or off vibrational entropy.  Add 3kB per added atom, remove that for removed atoms.  G = H-TS so interstiitals should decrease G while vacacnies add it
conditions.vib_ent_flag = 'On';
% conditions.vib_ent_flag = 'Off';

% Turn on or off T dependent band parameters
conditions.T_dep_bands_flag = 'On';
% conditions.T_dep_bands_flag = 'Off';

% set how to handle chemical potentials.  element order is defined as Ga, O, Si, H, Fe ....
conditions.T_dep_mu_flag = 'On';
%conditions.T_dep_mu_flag = 'Off';

conditions.impurity_solubility_flag = 'On';
% conditions.impurity_solubility_flag = 'Off';

% set how to handle defects with fixed concentrtations - either one
% constant value for all T or different values vs T
% conditions.T_dep_fixed_defect_flag = 'On';
conditions.T_dep_fixed_defect_flag = 'Off';

% choose the method of optimization used for finding solutions when elements are constrained - case 3 and 4
% conditions.search_method_flag = 'grid_fminsearch';
% conditions.search_method_flag = 'particleswarm_pattern';
conditions.search_method_flag = 'particleswarm_pattern_simplex';

if strcmp(conditions.search_method_flag,'grid_fminsearch')
    %% set options for fminsearch for fixing elements concentrations
    conditions.fixed_element_fmin_MaxFunEvals = 1e6;
    conditions.fixed_element_fmin_MaxIter = 1e6;
    conditions.fixed_element_fmin_TolX = 1e-5;
    conditions.fixed_element_fmin_TolFun = 1e-5;

elseif strcmp(conditions.search_method_flag,'particleswarm_pattern') || strcmp(conditions.search_method_flag,'particleswarm_pattern_simplex')
    conditions.fixed_element_swarm_iterlimit_base = 20;
    conditions.fixed_element_swarm_stall_lim = 15;
    conditions.fixed_element_swarm1_max_size = 15000;   %15000 seems to be ok upper limit on a laptop, 20000 is slow
    conditions.fixed_element_swarm2_search_band_kB = 3;
    conditions.fixed_element_swarm2_min_size = 3000;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






% 
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












    %% set total pressure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set total pressure
    % conditions.P_tot = 1e-4;   % total pressure in atm
    conditions.P_tot = 1;   % total pressure in atm
    conditions.P_units = 'atm';


    %% set chemical potentials for host material and impuritities %%%%%%%%%%%%%
    % note: later below some elements or individual defects' concentrations can
    % be fixed.  This overrides setting chemical potentials up here so it's OK
    % to have entries here for elements fixed later.

    %%% different options for setting the host material's chemical potentials according to equilibrium equations or other
    conditions.mu_conditions_flag = 'pO2-variable';   % equilibrium with pO2 set to arbitrary value.  you must also set pO2.
    % conditions.mu_conditions_flag = 'O2-1atm';   % pure O2 gas at 1 atm
    % conditions.mu_conditions_flag = 'pGa_vap-variable';   % equilibrium with pGa_vapor set to arbitrary value.  you must also set pGa_vap.
    % conditions.mu_conditions_flag = 'Ga_equilib_vap_pressure';    % set pGa_vap to equilibrium vapor pressure of Ga over Ga
    % conditions.mu_conditions_flag = 'Ga_ls-rich';  % use this to mean Ga2O3 in direct contact with pure Ga liq or solid 
    % conditions.mu_conditions_flag = 'pGa2O_vap-variable';   
    % conditions.mu_conditions_flag = 'pGaO_vap-variable';   
    
    
    %% conditions.mu_conditions_flag = 'manual' - add this sometime to account
    %% for things like growth or decomposition - excess chemical potential -
    %% here one or the other elements' mu would be greater than or less than the
    %% amount needed for equilibrium

    % set pO2 - if using this, need to also set conditions.mu_conditions_flag to "pO2-variable"
    % conditions.pO2 = 1;
    % conditions.pO2 = 0.2;
    conditions.pO2 = 0.02;
    % conditions.pO2 = 1e-4;
    % conditions.pO2 = 1e-6;
    % conditions.pO2 = 1e-8;
    % conditions.pO2 = 1e-10;

    % set pGa_vapor - if using this, need to also set conditions.mu_conditions_flag to "pGa_vap-variable"
    % conditions.pGa_vap = 1;
    % conditions.pGa_vap = 1e-4;
    % conditions.pGa_vap = 1e-8;
    % conditions.pGa_vap = 1e-10;

    % specify partial pressures of Ga2O or GaO
    % conditions.pGa2O_vap = 1e-4;
    % conditions.pGaO_vap = 1e-4;

    % set chemical potentials for some elements as constants
    % elements are in order Ga O Si H Fe Sn
    % Ga and O chemical potentials handled below - set only those for
    % impurities here.  If any elements are fixed, these values are ignored
    % (fixed elements supeercede these fixed mu values)

    % mu_Si = -20;
    % mu_H = -5;    % can link this to the pO2 through the H2O equilibrium
    % mu_Fe = -20;
    % mu_Sn = -4.5;  % works for pO2=1e-4 at high T to give 2e18 doping from SnGaII
    % mu_Sn = -5.3; %works for pO2=1 atm at 1000 C

    % set the mu values held constant here.  Any that are given fixed concentrations below will override these values.
    %%%%%% TO DO idea - actually set the unwated defects concentrations to zero
    %%%%%% - or use mu= -Inf?
    %%%%%%
    mu_Si = -30;
    % mu_H = -1.5;
    mu_H = -30;
    % mu_Fe = -5.5;   %This gives about 1e16 to 1e17
    mu_Fe = -30;
    mu_Sn = -30;
    % mu_Sn = -5; % need to change this value for every temperature to get the constant Sn concentration
    mu_Cr = -30;
    % mu_Cr = -8;  %will give about 1e16
    mu_Ti = -30;
    mu_Ir = -30;
    mu_Mg = -30;
    mu_Ca = -30;
    mu_Zn = -30;
    mu_Co = -30;
    mu_Zr = -30;
    mu_Hf = -30;
    mu_Ta = -30;
    mu_Ge = -30;
    impurity_mu_holder = [mu_Si mu_H mu_Fe mu_Sn mu_Cr mu_Ti mu_Ir mu_Mg mu_Ca mu_Zn mu_Co mu_Zr mu_Hf mu_Ta mu_Ge];    % this holds the mu for all the impurity elements - matrix Ga and O set below and added to this



    %%%%%%%%% can add some of these into defects database
    defects.mu_names = ["\mu_{Ga}", "\mu_{O}", "\mu_{Si}", "\mu_{H}", "\mu_{Fe}", "\mu_{Sn}", "\mu_{Cr}", "\mu_{Ti}", "\mu_{Ir}", "\mu_{Mg}", "\mu_{Ca}", "\mu_{Zn}", "\mu_{Co}", "\mu_{Zr}", "\mu_{Hf}", "\mu_{Ta}", "\mu_{Ge}"];
    defects.mu_names_with_units = ["\mu_{Ga} (eV)", "\mu_{O} (eV)", "\mu_{Si} (eV)", "\mu_{H} (eV)", "\mu_{Fe} (eV)", "\mu_{Sn} (eV)", "\mu_{Cr} (eV)", "\mu_{Ti} (eV)", "\mu_{Ir} (eV)", "\mu_{Mg} (eV)", "\mu_{Ca} (eV)", "\mu_{Zn} (eV)", "\mu_{Co} (eV)", "\mu_{Zr} (eV)", "\mu_{Hf} (eV)", "\mu_{Ta} (eV)", "\mu_{Ge} (eV)"];



    clear mu_Si mu_H mu_Fe mu_Sn mu_Cr mu_Ti mu_Ir mu_Mg mu_Ca mu_Zn mu_Co mu_Zr mu_Hf mu_Ta mu_Ge

    % conditions.num_elements = 17;

    % if conditions.num_elements is not explicitly set, set it to the number of
    % columns in defects.cs_dm
    if ~isfield(conditions,'num_elements') || sum(isempty(conditions.num_elements))
        conditions.num_elements = size(defects.cs_dm,2);
    end

    % TO DO - the order of the elements is idiosyncratic right now - we should
    % just make it the order of the periodic table for generality, so just make
    % cs_dm 118 columns wide etc.

    %%% select the element you want to vary mu for at fixed T for Task 5
    % conditions.Task5_element_num_to_vary_mu = 1; % for Ga
    % conditions.Task5_element_num_to_vary_mu = 2; % for O
    % conditions.Task5_element_num_to_vary_mu = 3; % for Si
    % conditions.Task5_element_num_to_vary_mu = 4; % for H
    % conditions.Task5_element_num_to_vary_mu = 5; % for Fe
    % conditions.Task5_element_num_to_vary_mu = 6; % for Sn

    %%%%%%%% end setting pressures and chemical potnetials  %%%%%%%%%%%%%%%%%%%





    %% set temperatures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % vector of T values at which to calculate for equilibrium
    % start at the highest T you want and then work down in order to fully
    % exploit the adaptive guessing.  Also, for fixed elements the major time
    % is spent finding the first guess, the optimization is fast by comparison once you fix 2 or more elements.
    % So it's worth doing a fine T grid for long calcuations so you dont have
    % to ever do it over.

    % conditions.T_equilibrium = 1500:-100:700;  % row vector of T in K at which to compute the equilirium
    conditions.T_equilibrium = 2100:-25:500;
    % conditions.T_equilibrium = 2100:-100:2000;
    % conditions.T_equilibrium = 10:10:300;
    % conditions.T_equilibrium = 2100:-50:1400;
    % conditions.T_equilibrium = 500:50:2100;
    % conditions.T_equilibrium = 500:100:2100;
    % conditions.T_equilibrium = 300:100:2100;
    % conditions.T_equilibrium = 1000 + 273;   % use this for tasks like doping or mu sweeps done at one temperature
    % conditions.T_equilibrium = 600 + 273;
    % conditions.T_equilibrium = 1600 + 273;
    % conditions.T_equilibrium = 1200 + 273;
    % conditions.T_equilibrium = 1100 + 273;
    % conditions.T_equilibrium = 1600 + 273;
    conditions.kBT_equilibrium = conditions.kB * conditions.T_equilibrium;


    % set final T for full quenching from the T_equilibrium values above
    conditions.T_fullquench = 300;   % scalar temperature at which to compute the defect concentrations after quenching from T_equilibrium
    conditions.kBT_fullquench = conditions.kB * conditions.T_fullquench;
    %%% end setting temperatures  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





    %% set unidentified shallow doping (excess charge added to charge balance)
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





    %% set any defects with fixed concentrations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% T independent case - one value for all temperatures
    % first initialize a COLUMNvector with the right number of entries for the
    % current defect variable.  Then you can go in and st individual ones
    conditions.fixed_defects = zeros(defects.num_defects,1);
    conditions.fixed_defects_concentrations = zeros(defects.num_defects,1);

    %
    % conditions.fixed_defects(169)= 1;
    % conditions.fixed_defects_concentrations(169) = 8.9e16;   % Si_GaI


     %% Sn doped Kuramata 2016
        conditions.fixed_defects(222)= 1;
        conditions.fixed_defects_concentrations(222) = 3.6e17;   % Ir_GaII
        conditions.fixed_defects(193)= 1;
        conditions.fixed_defects_concentrations(193) = 4.3e16;   % Cr_GaII
        conditions.fixed_defects(199)= 1;
        conditions.fixed_defects_concentrations(199) = 2.5e15;   % Ti_GaII
        conditions.fixed_defects(247)= 1;
        conditions.fixed_defects_concentrations(247) = 6.3e15;   % Zr_GaII
        conditions.fixed_defects(238)= 1;
        conditions.fixed_defects_concentrations(238) = 1.3e16;   % Ca_GaII
        conditions.fixed_defects(236)= 1;
        conditions.fixed_defects_concentrations(236) = 1.4e16;   % Mg_GaII













    % 
    % 
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
    % 
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
    % 
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
    % 
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






    %% Set elements with fixed concentrations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if an element is fixed, it means the chemical potential for that element
    % is floated (rather than fixed) and optimized to get the desired ttoal
    % number over all the defects containing that elemen
    % list any fixed elements by putting a 1 in its position in the vector conditions.fixed_elements

    % element order for ref
    % 1=Ga 2=O 3=Si 4=H 5=Fe 6=Sn 7=Cr 8=Ti 9=Ir 10=Mg 11=Ca 12=Zn 13=Co 14=Zr 15=Hf 16=Ta 17=Ge
    % conditions.fixed_elements = zeros(1,conditions.num_elements);
    % conditions.fixed_elements_concentrations = zeros(1,conditions.num_elements);

    % first initialize a vector with the right number of entries for the
    % current defect variable.  Then you can go in and st individual ones
    conditions.fixed_elements = zeros(1,defects.numelements);
    conditions.fixed_elements_concentrations = zeros(1,defects.numelements);


    % conditions.fixed_elements  = [0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];  % fix Si
    % conditions.fixed_elements  = [0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0];  % fix Sn
    % conditions.fixed_elements  = [0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0];  % fix Fe
    % conditions.fixed_elements_concentrations = [0 0 1e16 1e15 1e17 2e18 0 0 0 0 0 0 0 0 0 0 0];  % enter the conc of any fixed elements.  There can be nonzero values in here if the element is not fixed
    % 
    % 
    % %     % Sn KURAMATA 2016
    % 
    %     conditions.fixed_elements(6) = 1;    % fix Sn
    %     conditions.fixed_elements_concentrations(6) = 4.5e18;
    % 
    % 
    % 




    % 
    % 
    % if index==1 || index==5
    % 
    %     % Sn KURAMATA 2016
    %     conditions.fixed_elements(3) = 1;    % fix Si
    %     conditions.fixed_elements_concentrations(3) = 2.0e17;
    %     conditions.fixed_elements(6) = 1;    % fix Sn
    %     conditions.fixed_elements_concentrations(6) = 4.5e18;
    %     conditions.fixed_elements(5) = 1;    % fix Fe
    %     conditions.fixed_elements_concentrations(5) = 6.4e16;
    % 
    % elseif index==2 || index==6
    % 
    %     % Mg
    %     conditions.fixed_elements(3) = 1;    % fix Si
    %     conditions.fixed_elements_concentrations(3) = 2.3e17;
    %     conditions.fixed_elements(6) = 1;    % fix Sn
    %     conditions.fixed_elements_concentrations(6) = 6.3e15;
    %     conditions.fixed_elements(5) = 1;    % fix Fe
    %     conditions.fixed_elements_concentrations(5) = 3.0e16;
    % 
    % elseif index==3 || index==7
    % 
    %     % Fe
    %     conditions.fixed_elements(3) = 1;    % fix Si
    %     conditions.fixed_elements_concentrations(3) = 2.3e17;
    %     conditions.fixed_elements(6) = 1;    % fix Sn
    %     conditions.fixed_elements_concentrations(6) = 6.5e15;
    %     conditions.fixed_elements(5) = 1;    % fix Fe
    %     conditions.fixed_elements_concentrations(5) = 2.5e18;
    % 
    % elseif index==4 || index==8
    %     % UID
    %     conditions.fixed_elements(3) = 1;    % fix Si
    %     conditions.fixed_elements_concentrations(3) = 2.3e17;
    %     conditions.fixed_elements(6) = 1;    % fix Sn
    %     conditions.fixed_elements_concentrations(6) = 6.5e15;
    %     conditions.fixed_elements(5) = 1;    % fix Fe
    %     conditions.fixed_elements_concentrations(5) = 3.0e16;
    % 
    % end
    % 






    % conditions.fixed_elements(4) = 1;    % fix H at 1e17
    % conditions.fixed_elements_concentrations(4) = 1e17;
    %
    % conditions.fixed_elements(3) = 1;    % fix Si at 1e17
    % conditions.fixed_elements_concentrations(3) = 5e16;
    %
    % conditions.fixed_elements(3) = 1;    % fix Si at 8.9e16
    % conditions.fixed_elements_concentrations(3) = 8.9e16;
    %
    % conditions.fixed_elements(4) = 1;    % fix H at 1e17
    % conditions.fixed_elements_concentrations(4) = 1e17;

    % conditions.fixed_elements(5) = 1;    % fix Fe
    % conditions.fixed_elements_concentrations(5) = 2.8e16;

    % conditions.fixed_elements(6) = 1;    % fix Sn at 2e18
    % conditions.fixed_elements_concentrations(6) = 2e18;
    %
    % conditions.fixed_elements(7) = 1;    % fix Cr
    % conditions.fixed_elements_concentrations(7) = 1.9e16;

    % conditions.fixed_elements(8) = 1;    % fix Ti at 1e16
    % conditions.fixed_elements_concentrations(8) = 1e16;

    % conditions.fixed_elements(9) = 1;    % fix Ir at 1.7e17
    % conditions.fixed_elements_concentrations(9) = 1.7e17;

    % set the ranges for each mu to search over when elements are fixed.
    % since a direc grid search is used to generate the first guess, there are
    % HUGE time savings if you can narrow down this range manually first using
    % a series of calcs with constant mu values to estimate mu for each element
    % you want to freeze.  Or you can brute force it and just wait for ever...


    % there need to be values here for the ranges

    % lo_mu_Ga = 0;
    % hi_mu_Ga = 0;
    %
    % lo_mu_O = 0;
    % hi_mu_O = 0;
    % 
    % lo_mu_Si = -10;
    % hi_mu_Si = -2;
    % %
    % % lo_mu_H = -10;
    % % hi_mu_H = -2; % can link this to the pO2 through the H2O equilibrium
    % 
    % % 1e16 to 1e17 Fe is about -6 to -5 eV with 2e18 Sn doping
    % lo_mu_Fe = -10;
    % hi_mu_Fe = -1;
    % 
    % % 2e18 Sn is about -5.5 eV at 2100K  to -5 eV at 500 K
    % lo_mu_Sn = -10;
    % hi_mu_Sn = -2;

    % % 1e16 Cr is about -8 to -9 for 500-2100 K
    % lo_mu_Cr = -10;
    % hi_mu_Cr = -4;
    %
    % % 1e16 Ti
    % lo_mu_Ti = -10;
    % hi_mu_Ti = -5;
    %
    % % 1e17 Ir
    % lo_mu_Ir = -3;
    % hi_mu_Ir = -1.5;



    % construct the ranges over which to search for mu for frozen elements.  convention for now is that this matrix should be the size for all the elements, not just the fixed ones we specify ranges for all the elements.  Later can change this to only specify ranges for the ones that are frozen (not specify values for those left open).
    conditions.fixed_elements_mu_ranges = zeros(conditions.num_elements,2);   % create  matrix of right size for ALL the elements - only will pay attention to those which are fixed




    %% improve this by doing these steps automatically if it's detected that the ith element is frozen.  Then do a logic check to make sure upper and lower bounds are given and if none given, choose a default
    % Ga and O are 1 and 2
    % conditions.fixed_elements_mu_ranges(3,:) = [lo_mu_Si hi_mu_Si];
    % conditions.fixed_elements_mu_ranges(4,:) = [lo_mu_H hi_mu_H];
    % conditions.fixed_elements_mu_ranges(5,:) = [lo_mu_Fe hi_mu_Fe];
    % conditions.fixed_elements_mu_ranges(6,:) = [lo_mu_Sn hi_mu_Sn];
    % conditions.fixed_elements_mu_ranges(7,:) = [lo_mu_Cr hi_mu_Cr];
    % conditions.fixed_elements_mu_ranges(8,:) = [lo_mu_Ti hi_mu_Ti];
    % conditions.fixed_elements_mu_ranges(9,:) = [lo_mu_Ir hi_mu_Ir];


    %% add in ability to import T-dependnet Ef and mu's from a prior simulation.  Let's say you fixed 2 elements, got a good solution, and next want to freeze a 3rd one.  The best place to start is that last solution then ramp up the fixed concentration on the new element from 0 to it's target value.
    %% actually could make it automatically do this for multi-element problems - start with the largest target concentration element, then add the next a perturbing one especially if it forms complexes with that first one, then go to the next element..





    % default is no elements fixed - catch if nothing entred.
    if ~isfield(conditions,'fixed_elements') || sum(isempty(conditions.fixed_elements))
        conditions.fixed_elements = zeros(defects.numelements,1);
    end

    if ~isfield(conditions,'fixed_elements_concentrations') || sum(isempty(conditions.fixed_elements_concentrations))
        conditions.fixed_elements_concentrations = zeros(defects.numelements,1);
    end

    % if something is not fixed but by accident it has a target conc set, jsut
    % set its target concentration to zero also, using the fixed defects vector
    % as a mask
    conditions.fixed_elements_concentrations = conditions.fixed_elements_concentrations.*conditions.fixed_elements;

    % find number of fixed elements and their indices
    conditions.indices_of_fixed_elements = find(conditions.fixed_elements);
    conditions.num_fixed_elements = numel(conditions.indices_of_fixed_elements);


    if conditions.num_fixed_elements~=0
        for i=1:conditions.num_fixed_elements
            msg = strcat('WARNING: Fixed concentration of element used for...',defects.elementnames(conditions.indices_of_fixed_elements(i)),'...Be sure you meant to do this.');
            disp(msg)
        end
    elseif conditions.num_fixed_defects==0
        disp('No Eelements fixed.  Proceeding')
    end

    clear i msg lo_mu_Ga hi_mu_Ga lo_mu_O hi_mu_O lo_mu_Si hi_mu_Si lo_mu_H hi_mu_H lo_mu_Fe hi_mu_Fe lo_mu_Sn hi_mu_Sn lo_mu_Cr hi_mu_Cr lo_mu_Ti hi_mu_Ti lo_mu_Ir hi_mu_Ir max_mu_range start_kBT fine_factor
    %%% end section on fixed elements %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




    %% Set Semiconductor Bandstructure and Lattice sites %%%%%%%%%%%%%%%%%%%%%%
    %%% Ga2O3 %%
    conditions.TRef = 300;
    conditions.EgRef = 4.8;   % Varshini equation or other Eg(T) parameterization should be consistent with this
    conditions.EvRef = 0;
    conditions.EcRef = conditions.EgRef;
    conditions.NcRef = 1e19;
    conditions.NvRef = 5e20;
    conditions.EcT_fraction = 0.75;     % what fraction of Eg change vs T happens in the CB?
    % conditions.EcT_fraction = 0.50;     % what fraction of Eg change vs T happens in the CB?
    conditions.EvT_fraction = 1-conditions.EcT_fraction;
    conditions.Eg0 = 4.9;   % bandgap at 0K in eV.  Using Adrian's number from STEM EELS at PSU - 4.8 eV at 300 K.  So about 4.9 eV at 0K using
    conditions.varshini_a = 0.00105248259206626;
    conditions.varshini_b = 676.975385507958;

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




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% handle setting matrix chemical potentials based on inputs above this line
    % note that the mu_values for the matrix elements are T dependent so this
    % ends up creatiung a matrix not just a vector
    [mu_Ga, mu_O] = Ga_O_chem_potentials_from_conditions(conditions);
    % % [mu_Ga, mu_O] = Ga_O_chem_potentials_from_conditions(conditions,'yes');
    temp_mu = ones(size(conditions.T_equilibrium,2),1)*impurity_mu_holder;  % this has to be here or else horzcat throws an error - creates a matrix with constant mu values for each element at each tempetature by outer product of temp vector and mu values set for impurities
    conditions.muT_equilibrium = [mu_Ga mu_O temp_mu];  % commits all the chem potentials to the conditions variable
    if size(conditions.muT_equilibrium,2)~=conditions.num_elements
        error('number of columns in mu array not qual to number of elements declared')
    end
    clear mu_Ga mu_O temp_mu impurity_mu_holder
    %%% end handling matrix chemical potentials %%%%%%%%%%%%

    if strcmp(conditions.impurity_solubility_flag,'On')
        % [conditions.muT_equilibrium(:,3), conditions.muT_equilibrium(:,4), conditions.muT_equilibrium(:,5), conditions.muT_equilibrium(:,6), ~, ~, conditions.muT_equilibrium(:,9), conditions.muT_equilibrium(:,10), conditions.muT_equilibrium(:,11)] = mu_limit_imposed_by_pO2(conditions.T_equilibrium, conditions.P_tot, conditions.P_units, conditions.muT_equilibrium(:,2));
        % [conditions.muT_equilibrium(:,3), ~, conditions.muT_equilibrium(:,5), conditions.muT_equilibrium(:,6), ~, ~, ~, conditions.muT_equilibrium(:,10), conditions.muT_equilibrium(:,11), ~, ~, ~, ~, ~, ~] = mu_limit_imposed_by_pO2(conditions.T_equilibrium, conditions.P_tot, conditions.P_units, conditions.muT_equilibrium(:,2));
        [conditions.muT_equilibrium(:,3), ~, conditions.muT_equilibrium(:,5), conditions.muT_equilibrium(:,6), ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = mu_limit_imposed_by_pO2(conditions.T_equilibrium, conditions.P_tot, conditions.P_units, conditions.muT_equilibrium(:,2));
        
        % [conditions.muT_equilibrium(:,3), ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = mu_limit_imposed_by_pO2(conditions.T_equilibrium, conditions.P_tot, conditions.P_units, conditions.muT_equilibrium(:,2));
        % [~, conditions.muT_equilibrium(:,4), ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = mu_limit_imposed_by_pO2(conditions.T_equilibrium, conditions.P_tot, conditions.P_units, conditions.muT_equilibrium(:,2));
        % [~, ~, conditions.muT_equilibrium(:,5), ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = mu_limit_imposed_by_pO2(conditions.T_equilibrium, conditions.P_tot, conditions.P_units, conditions.muT_equilibrium(:,2));
        % [~, ~, ~, conditions.muT_equilibrium(:,6), ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = mu_limit_imposed_by_pO2(conditions.T_equilibrium, conditions.P_tot, conditions.P_units, conditions.muT_equilibrium(:,2));
        % [~, ~, ~, ~,conditions.muT_equilibrium(:,7),  ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = mu_limit_imposed_by_pO2(conditions.T_equilibrium, conditions.P_tot, conditions.P_units, conditions.muT_equilibrium(:,2));
        % [~, ~, ~, ~, ~, conditions.muT_equilibrium(:,8), ~, ~, ~, ~, ~, ~, ~, ~, ~] = mu_limit_imposed_by_pO2(conditions.T_equilibrium, conditions.P_tot, conditions.P_units, conditions.muT_equilibrium(:,2));
        % [~, ~, ~, ~, ~, ~, conditions.muT_equilibrium(:,9), ~, ~, ~, ~, ~, ~, ~, ~] = mu_limit_imposed_by_pO2(conditions.T_equilibrium, conditions.P_tot, conditions.P_units, conditions.muT_equilibrium(:,2));
        % [~, ~, ~, ~, ~, ~, ~,conditions.muT_equilibrium(:,10), ~, ~, ~, ~, ~, ~, ~] = mu_limit_imposed_by_pO2(conditions.T_equilibrium, conditions.P_tot, conditions.P_units, conditions.muT_equilibrium(:,2));
        % [~, ~, ~, ~, ~, ~, ~, ~,conditions.muT_equilibrium(:,11), ~, ~, ~, ~, ~, ~] = mu_limit_imposed_by_pO2(conditions.T_equilibrium, conditions.P_tot, conditions.P_units, conditions.muT_equilibrium(:,2));
        % [~, ~, ~, ~, ~, ~, ~, ~, ~,conditions.muT_equilibrium(:,12), ~, ~, ~, ~, ~] = mu_limit_imposed_by_pO2(conditions.T_equilibrium, conditions.P_tot, conditions.P_units, conditions.muT_equilibrium(:,2));
        % [~, ~, ~, ~, ~, ~, ~, ~, ~, ~,conditions.muT_equilibrium(:,13), ~, ~, ~, ~] = mu_limit_imposed_by_pO2(conditions.T_equilibrium, conditions.P_tot, conditions.P_units, conditions.muT_equilibrium(:,2));
        % [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~,conditions.muT_equilibrium(:,14), ~, ~, ~] = mu_limit_imposed_by_pO2(conditions.T_equilibrium, conditions.P_tot, conditions.P_units, conditions.muT_equilibrium(:,2));
        % [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~,conditions.muT_equilibrium(:,15), ~, ~] = mu_limit_imposed_by_pO2(conditions.T_equilibrium, conditions.P_tot, conditions.P_units, conditions.muT_equilibrium(:,2));
        % [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, conditions.muT_equilibrium(:,16), ~] = mu_limit_imposed_by_pO2(conditions.T_equilibrium, conditions.P_tot, conditions.P_units, conditions.muT_equilibrium(:,2));
        % [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~,conditions.muT_equilibrium(:,17)] = mu_limit_imposed_by_pO2(conditions.T_equilibrium, conditions.P_tot, conditions.P_units, conditions.muT_equilibrium(:,2));

        % % % Si=3  H=4  Fe=5  Sn=6  Cr=7  Ti=8  Ir=9  Mg=10  Ca=11  Zn=12  Co=13  Zr=14  Hf=15  Ta=16  Ge=17


    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if sum(sum(isinf(conditions.muT_equilibrium)))
        error('Infinite mu detected')
    end


    %% save the problem definition in the conditions variable into a mat file for later inspection
    save(strcat(save_pname,conditions_save_fname,'.mat'),'conditions')
    disp('Saved problem specificaiton conditions to .mat file')
    clear conditions_save_fname
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
    [equilib_dark_sol] = defect_equilibrium_dark_with_frozen(conditions, defects);
    toc

    
    %%%%%%%%%%% save the outputs
    all_defects_out = cat(2,equilib_dark_sol.T_equilibrium, equilib_dark_sol.charge_bal_err, equilib_dark_sol.element_bal_err, equilib_dark_sol.tot_bal_err, equilib_dark_sol.mu, conditions.EvT_equilibrium', conditions.EcT_equilibrium', equilib_dark_sol.EFn, equilib_dark_sol.EFp, equilib_dark_sol.Nd, equilib_dark_sol.Na, equilib_dark_sol.n, equilib_dark_sol.p, equilib_dark_sol.defects);
    all_chargestates_out = cat(2,equilib_dark_sol.T_equilibrium, equilib_dark_sol.charge_bal_err, equilib_dark_sol.element_bal_err, equilib_dark_sol.tot_bal_err, equilib_dark_sol.mu, conditions.EvT_equilibrium', conditions.EcT_equilibrium', equilib_dark_sol.EFn, equilib_dark_sol.EFp, equilib_dark_sol.Nd, equilib_dark_sol.Na, equilib_dark_sol.n, equilib_dark_sol.p, equilib_dark_sol.chargestates);
    all_defect_headers = ['T_equilibrium (K)' 'charge_bal_err' 'element_bal_err' 'tot_bal_err' defects.mu_names_with_units 'Ev(T) (eV)' 'Ec(T) (eV)' 'Efn(eV)' 'Efp (eV)' 'Nd (#/cm3)' 'Na (#/cm3)' 'n (#/cm3)' 'p (#/cm3)' defects.defect_names'];
    all_chargestate_headers = ['T_equilibrium (K)' 'charge_bal_err' 'element_bal_err' 'tot_bal_err' defects.mu_names_with_units 'Ev(T) (eV)' 'Ec(T) (eV)' 'Efn(eV)' 'Efp (eV)' 'Nd (#/cm3)' 'Na (#/cm3)' 'n (#/cm3)' 'p (#/cm3)' defects.chargestate_names'];
    all_defects_cell = [cellstr(all_defect_headers); num2cell(all_defects_out)];
    all_chargestates_cell = [cellstr(all_chargestate_headers); num2cell(all_chargestates_out)];

    


    % throw out all the ones with insignificant concentrations and save
    sig_defects_index = find(max(equilib_dark_sol.defects,[],1)>=conditions.save_min);
    sig_chargestates_index = find(max(equilib_dark_sol.chargestates,[],1)>=conditions.save_min);
    sig_defects_out = cat(2,equilib_dark_sol.T_equilibrium, equilib_dark_sol.charge_bal_err, equilib_dark_sol.element_bal_err, equilib_dark_sol.tot_bal_err, equilib_dark_sol.mu, conditions.EvT_equilibrium', conditions.EcT_equilibrium', equilib_dark_sol.EFn, equilib_dark_sol.EFp, equilib_dark_sol.Nd, equilib_dark_sol.Na, equilib_dark_sol.n, equilib_dark_sol.p, equilib_dark_sol.defects(:,sig_defects_index));
    sig_chargestates_out = cat(2,equilib_dark_sol.T_equilibrium, equilib_dark_sol.charge_bal_err, equilib_dark_sol.element_bal_err, equilib_dark_sol.tot_bal_err, equilib_dark_sol.mu, conditions.EvT_equilibrium', conditions.EcT_equilibrium', equilib_dark_sol.EFn, equilib_dark_sol.EFp, equilib_dark_sol.Nd, equilib_dark_sol.Na, equilib_dark_sol.n, equilib_dark_sol.p, equilib_dark_sol.chargestates(:,sig_chargestates_index));
    sig_defect_headers = ['T_equilibrium (K)' 'charge_bal_err' 'element_bal_err' 'tot_bal_err' defects.mu_names_with_units 'Ev(T) (eV)' 'Ec(T) (eV)' 'Efn(eV)' 'Efp (eV)' 'Nd (#/cm3)' 'Na (#/cm3)' 'n (#/cm3)' 'p (#/cm3)' defects.defect_names(sig_defects_index)'];
    sig_chargestate_headers = ['T_equilibrium (K)' 'charge_bal_err' 'element_bal_err' 'tot_bal_err' defects.mu_names_with_units 'Ev(T) (eV)' 'Ec(T) (eV)' 'Efn(eV)' 'Efp (eV)' 'Nd (#/cm3)' 'Na (#/cm3)' 'n (#/cm3)' 'p (#/cm3)' defects.chargestate_names(sig_chargestates_index)'];
    sig_defects_cell = [cellstr(sig_defect_headers); num2cell(sig_defects_out)];
    sig_chargestates_cell = [cellstr(sig_chargestate_headers); num2cell(sig_chargestates_out)];

    % write output tables to a tab delim text file
    % writematrix(equilib_dark_sol_out,strcat(save_pname,equilib_save_fname),'Delimiter','tab','WriteMode','overwrite')
    writecell(all_defects_cell,strcat(save_pname,equilib_save_fname,'_all_defects'),'Delimiter','tab','WriteMode','overwrite');  % write results with column headers
    writecell(all_chargestates_cell,strcat(save_pname,equilib_save_fname,'_all_chargestates'),'Delimiter','tab','WriteMode','overwrite');
    writecell(sig_defects_cell,strcat(save_pname,equilib_save_fname,'_sig_defects'),'Delimiter','tab','WriteMode','overwrite');
    writecell(sig_chargestates_cell,strcat(save_pname,equilib_save_fname,'_sig_chargestates'),'Delimiter','tab','WriteMode','overwrite');
    disp('Wrote equilibrium dark solution as four text files')

    % write the whole solution structure to a mat file
    save(strcat(save_pname,equilib_save_fname))
    disp('Wrote equilibrium dark solution to .mat file')

    % clean up memory
    clear equilib_save_fname all_defects_out all_chargestates_out all_defects_headers all_chargestates_headers all_defects_cell all_chargstates_cell sig_defects_index sig_chargestates_index sig_defects_out sig_chargestates_out sig_defects_headers sig_chargestates_headers sig_defects_cell sig_chargestates_cell
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    % %% Task Trange 2: quenching from the annealing temperatures from example 1  %%%%%%%%%%%%%%%
    % 
    % %%% note this assumes equilib_dark_sol is still in memory - could also
    % %%% add a line here to load a prior solution from a mat file.
    % 
    % tic
    % [fullquench_dark_sol] = defect_fullquench_dark(equilib_dark_sol, conditions, defects);
    % toc
    % 
    % % these variables don't exist for quenching so set them to zero to keep the
    % % file format the same with the full equilibrium calc
    % fullquench_dark_sol.element_bal_err = zeros(size(fullquench_dark_sol.charge_bal_err));
    % fullquench_dark_sol.tot_bal_err = zeros(size(fullquench_dark_sol.charge_bal_err));
    % fullquench_dark_sol.mu = zeros(size(equilib_dark_sol.mu));
    % EvT_fullquench_vec = conditions.EvT_fullquench*ones(size(fullquench_dark_sol.charge_bal_err));
    % EcT_fullquench_vec = conditions.EcT_fullquench*ones(size(fullquench_dark_sol.charge_bal_err));
    % 
    % %%%%%%%%%%% save the outputs for quenching.  Yes some of these are the
    % %%%%%%%%%%% same as computed for the equilibrium solution but for sake
    % %%%%%%%%%%% of debugging/proofreading were jsut doing it here again
    % all_defects_out = cat(2,fullquench_dark_sol.T_equilibrium, fullquench_dark_sol.charge_bal_err, fullquench_dark_sol.element_bal_err, fullquench_dark_sol.tot_bal_err, fullquench_dark_sol.mu, conditions.EvT_equilibrium', conditions.EcT_equilibrium', fullquench_dark_sol.EFn, fullquench_dark_sol.EFp, fullquench_dark_sol.Nd, fullquench_dark_sol.Na, fullquench_dark_sol.n, fullquench_dark_sol.p, fullquench_dark_sol.defects);
    % all_chargestates_out = cat(2,fullquench_dark_sol.T_equilibrium, fullquench_dark_sol.charge_bal_err, fullquench_dark_sol.element_bal_err, fullquench_dark_sol.tot_bal_err, fullquench_dark_sol.mu, conditions.EvT_equilibrium', conditions.EcT_equilibrium', fullquench_dark_sol.EFn, fullquench_dark_sol.EFp, fullquench_dark_sol.Nd, fullquench_dark_sol.Na, fullquench_dark_sol.n, fullquench_dark_sol.p, fullquench_dark_sol.chargestates);
    % all_defect_headers = ['T_equilibrium (K)' 'charge_bal_err' 'element_bal_err' 'tot_bal_err' defects.mu_names_with_units 'Ev(T) (eV)' 'Ec(T) (eV)' 'Efn(eV)' 'Efp (eV)' 'Nd (#/cm3)' 'Na (#/cm3)' 'n (#/cm3)' 'p (#/cm3)' defects.defect_names'];
    % all_chargestate_headers = ['T_equilibrium (K)' 'charge_bal_err' 'element_bal_err' 'tot_bal_err' defects.mu_names_with_units 'Ev(T) (eV)' 'Ec(T) (eV)' 'Efn(eV)' 'Efp (eV)' 'Nd (#/cm3)' 'Na (#/cm3)' 'n (#/cm3)' 'p (#/cm3)' defects.chargestate_names'];
    % all_defects_cell = [cellstr(all_defect_headers); num2cell(all_defects_out)];
    % all_chargestates_cell = [cellstr(all_chargestate_headers); num2cell(all_chargestates_out)];
    % 
    % % throw out all the ones with insignificant concentrations and save
    % sig_defects_index = find(max(fullquench_dark_sol.defects,[],1)>=conditions.save_min);
    % sig_chargestates_index = find(max(fullquench_dark_sol.chargestates,[],1)>=conditions.save_min);
    % sig_defects_out = cat(2,fullquench_dark_sol.T_equilibrium, fullquench_dark_sol.charge_bal_err, fullquench_dark_sol.element_bal_err, fullquench_dark_sol.tot_bal_err, fullquench_dark_sol.mu, conditions.EvT_equilibrium', conditions.EcT_equilibrium', fullquench_dark_sol.EFn, fullquench_dark_sol.EFp, fullquench_dark_sol.Nd, fullquench_dark_sol.Na, fullquench_dark_sol.n, fullquench_dark_sol.p, fullquench_dark_sol.defects(:,sig_defects_index));
    % sig_chargestates_out = cat(2,fullquench_dark_sol.T_equilibrium, fullquench_dark_sol.charge_bal_err, fullquench_dark_sol.element_bal_err, fullquench_dark_sol.tot_bal_err, fullquench_dark_sol.mu, conditions.EvT_equilibrium', conditions.EcT_equilibrium', fullquench_dark_sol.EFn, fullquench_dark_sol.EFp, fullquench_dark_sol.Nd, fullquench_dark_sol.Na, fullquench_dark_sol.n, fullquench_dark_sol.p, fullquench_dark_sol.chargestates(:,sig_chargestates_index));
    % sig_defects_headers = ['T_equilibrium (K)' 'charge_bal_err' 'element_bal_err' 'tot_bal_err' defects.mu_names_with_units 'Ev(T) (eV)' 'Ec(T) (eV)' 'Efn(eV)' 'Efp (eV)' 'Nd (#/cm3)' 'Na (#/cm3)' 'n (#/cm3)' 'p (#/cm3)' defects.defect_names(sig_defects_index)'];
    % sig_chargestate_headers = ['T_equilibrium (K)' 'charge_bal_err' 'element_bal_err' 'tot_bal_err' defects.mu_names_with_units 'Ev(T) (eV)' 'Ec(T) (eV)' 'Efn(eV)' 'Efp (eV)' 'Nd (#/cm3)' 'Na (#/cm3)' 'n (#/cm3)' 'p (#/cm3)' defects.chargestate_names(sig_chargestates_index)'];
    % sig_defects_cell = [cellstr(sig_defect_headers); num2cell(sig_defects_out)];
    % sig_chargestates_cell = [cellstr(sig_chargestate_headers); num2cell(sig_chargestates_out)];
    % 
    % % write to a tab delim text file
    % % writematrix(equilib_dark_sol_out,strcat(save_pname,equilib_save_fname),'Delimiter','tab','WriteMode','overwrite')
    % writecell(all_defects_cell,strcat(save_pname,fullquench_save_fname,'_all_defects'),'Delimiter','tab','WriteMode','overwrite');  % write results with column headers
    % writecell(all_chargestates_cell,strcat(save_pname,fullquench_save_fname,'_all_chargestates'),'Delimiter','tab','WriteMode','overwrite');
    % writecell(sig_defects_cell,strcat(save_pname,fullquench_save_fname,'_sig_defects'),'Delimiter','tab','WriteMode','overwrite');
    % writecell(sig_chargestates_cell,strcat(save_pname,fullquench_save_fname,'_sig_chargestates'),'Delimiter','tab','WriteMode','overwrite');
    % disp('Wrote equilibrium dark solution as four text files')
    % 
    % 
    % % write the whole solution structure to a mat file
    % save(strcat(save_pname,fullquench_save_fname))
    % disp('Wrote equilibrium dark solution to .mat file')
    % 
    % % clean up memory
    % clear fullquench_save_fname all_defects_out all_chargestates_out all_defects_headers all_chargestates_headers all_defects_cell all_chargstates_cell sig_defects_index sig_chargestates_index sig_defects_out sig_chargestates_out sig_defects_headers sig_chargestates_headers sig_defects_cell sig_chargestates_cell
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%% end quenched calc %%%%%%%%%%%%%%%%%






    %% Display the results


    %% plot total defect concentrations - these will be same for equilib and quenched except for n and p
    figure(2)
    f2line(1) = plot(conditions.T_equilibrium,log10(equilib_dark_sol.n),'o','DisplayName','n equilib','MarkerEdgeColor',"r",'MarkerFaceColor',"r");  % full circles for equilib, open for fullquench
    f2line(2) = plot(conditions.T_equilibrium,log10(equilib_dark_sol.p),'o','DisplayName','p equilib','MarkerEdgeColor',"b",'MarkerFaceColor',"b");
    f2line(3) = plot(conditions.T_equilibrium,log10(abs(equilib_dark_sol.Nd - equilib_dark_sol.Na)),'k-','DisplayName','|Nd-Na|');
    % f2line(4) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.n),'o','DisplayName','n quench','MarkerEdgeColor',"r");
    % f2line(5) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.p),'o','DisplayName','p quench','MarkerEdgeColor',"b");

  for i = 1:3
        f2line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f2line(i).DisplayName},size(f2line(i).XData)));
  end
legend('n equilib','p equilib','|Nd-Na|')

    % for i = 1:5
    %     f2line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f2line(i).DisplayName},size(f2line(i).XData)));
    % end
    % legend('n equilib','p equilib','|Nd-Na|', 'n quench', 'p quench')
    
    for i=1:defects.num_defects
        if max(log10(equilib_dark_sol.defects(:,i))) >= conditions.plot_log10_min   % only plot defects that at some point rise above 1e10 /cm3 - this will cut down the legend size
            f2line(i+5) = plot(conditions.T_equilibrium,log10(equilib_dark_sol.defects(:,i)),'DisplayName', defects.defect_names(i));
            f2line(i+5).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f2line(i+5).DisplayName},size(f2line(i+5).XData)));
        end
    end
    legend show
    title('Concentrations at Equilibrium Temperature K')
    xlabel('Equilibrium T (K)')
    ylabel('Log_{10} Concentrations (#/cm3)')
    ylim([conditions.plot_log10_min conditions.plot_log10_max])
    % xlim([conditions.T_equilibrium(1) conditions.T_equilibrium(end)])
    datacursormode.Enable='on';   %turn on the datatip to give name of defect on mouse hover




    %% plot all chargestates for equilibrium
    figure(3)
    f3line(1) = plot(conditions.T_equilibrium,log10(equilib_dark_sol.n),'o','DisplayName','n','MarkerEdgeColor',"r",'MarkerFaceColor',"r");
    f3line(2) = plot(conditions.T_equilibrium,log10(equilib_dark_sol.p),'o','DisplayName','p','MarkerEdgeColor',"b",'MarkerFaceColor',"b");
    f3line(3) = plot(conditions.T_equilibrium,log10(abs(equilib_dark_sol.Nd - equilib_dark_sol.Na)),'k-','DisplayName','|Nd-Na|');
    for i = 1:3
        f3line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f3line(i).DisplayName},size(f3line(i).XData)));
    end
    legend('n','p','|Nd-Na|')
    for i=1:defects.num_chargestates
        if max(log10(equilib_dark_sol.chargestates(:,i))) >= conditions.plot_log10_min   % only plot defects that at some point rise above 1e10 /cm3 - this will cut down the legend size
            f3line(i+3) = plot(conditions.T_equilibrium,log10(equilib_dark_sol.chargestates(:,i)),'DisplayName', defects.chargestate_names(i));
            f3line(i+3).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f3line(i+3).DisplayName},size(f3line(i+3).XData)));
        end
    end
    legend show
    title('Concentrations at Equilibrium Temperature K')
    xlabel('Equilibrium T (K)')
    ylabel('Log_{10} Concentrations (#/cm3)')
    ylim([conditions.plot_log10_min conditions.plot_log10_max])
    % xlim([conditions.T_equilibrium(1) conditions.T_equilibrium(end)])
    datacursormode.Enable='on';   %turn on the datatip to give name of defect on mouse hover



    % % plot all chargestates for quenching - chargestates get changed for quenching
    % figure(4)
    % f4line(1) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.n),'o','DisplayName','n','MarkerEdgeColor',"r");
    % f4line(2) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.p),'o','DisplayName','p','MarkerEdgeColor',"b");
    % f4line(3) = plot(conditions.T_equilibrium,log10(abs(fullquench_dark_sol.Nd - fullquench_dark_sol.Na)),'k-','DisplayName','|Nd-Na|');
    % for i = 1:3
    %     f4line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f4line(i).DisplayName},size(f4line(i).XData)));
    % end
    % legend('n','p','|Nd-Na|')
    % for i=1:defects.num_chargestates
    %     if max(log10(fullquench_dark_sol.chargestates(:,i))) >= conditions.plot_log10_min   % only plot defects that at some point rise above 1e10 /cm3 - this will cut down the legend size
    %         f4line(i+3) = plot(conditions.T_equilibrium,log10(fullquench_dark_sol.chargestates(:,i)),'DisplayName', defects.chargestate_names(i));
    %         f4line(i+3).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f4line(i+3).DisplayName},size(f4line(i+3).XData)));
    %     end
    % end
    % legend show
    % title('Concentrations at Equilibrium Temperature K')
    % xlabel('Equilibrium T (K)')
    % ylabel('Log_{10} Concentrations (#/cm3)')
    % ylim([conditions.plot_log10_min conditions.plot_log10_max])
    % xlim([conditions.T_equilibrium(1) conditions.T_equilibrium(end)])
    % datacursormode.Enable='on';   %turn on the datatip to give name of defect on mouse hover



    %% Plot the  band diagram vs T Ec, Ev, EFn and EFp vs T
    figure(5)

    plot(conditions.T_equilibrium,equilib_dark_sol.EFp,'b-')   % note in dark Efn and Efp will lie over eachother
    plot(conditions.T_equilibrium,equilib_dark_sol.EFn,'r-')
    % plot(conditions.T_equilibrium,fullquench_dark_sol.EFp,'--b')
    % plot(conditions.T_equilibrium,fullquench_dark_sol.EFn,'--r')

    %% handle the T dependent band edges or constant ones
    % if strcmp(conditions.T_dep_bands_flag,'On')
    plot(conditions.T_equilibrium,conditions.EcT_equilibrium,'k-')
    plot(conditions.T_equilibrium,conditions.EvT_equilibrium,'k-')
    % elseif strcmp(conditions.T_dep_bands_flag,'Off')
    plot(conditions.T_equilibrium,conditions.EcRef*ones(size(conditions.T_equilibrium)),'k--')
    plot(conditions.T_equilibrium,conditions.EvRef*ones(size(conditions.T_equilibrium)),'k--')
    % else
    % error('conditions.T_dep_bands_flag must be on or off')
    % end

    title('Fermi Levels for Equilibrium and Quenching vs T')
    legend('EFp(T_{Equilib})','EFn(T_{Equilib})', 'EFp(T_{Quench})', 'EFn(T_{Quench})','Ec(T_{Equilib})','Ev(T_{Equilib})','Ec(T_{Quench})','Ev(T_{Quench})')
    xlabel('Equilibrium T (K)')
    ylabel('Energies (eV)')





    %% Plot the mu values for the solution
    figure(6)

    %% handle the T dependent band edges or constant ones
    % for i=1:defects.numelements
    %     plot(conditions.T_equilibrium, equilib_dark_sol.mu)
    % end

    f6line(1) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,1),'DisplayName','Ga');
    f6line(2) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,2),'DisplayName','O');
    f6line(3) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,3),'DisplayName','Si');
    f6line(4) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,4),'DisplayName','H');
    f6line(5) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,5),'DisplayName','Fe');
    f6line(6) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,6),'DisplayName','Sn');
    f6line(7) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,7),'DisplayName','Cr');
    f6line(8) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,8),'DisplayName','Ti');
    f6line(9) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,9),'DisplayName','Ir');
    f6line(10) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,10),'DisplayName','Mg');
    f6line(11) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,11),'DisplayName','Ca');
    f6line(12) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,12),'DisplayName','Zn');
    f6line(13) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,13),'DisplayName','Co');
    f6line(14) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,14),'DisplayName','Zr');
    f6line(15) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,15),'DisplayName','Hf');
    f6line(16) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,16),'DisplayName','Ta');
    f6line(17) = plot(conditions.T_equilibrium,equilib_dark_sol.mu(:,17),'DisplayName','Ge');

    for i = 1:17
        f6line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f6line(i).DisplayName},size(f6line(i).XData)));
    end
    title('Chemical Potentials for Elements at Equilibrium and Quenching vs T')
    legend('\mu_{Ga}', '\mu_{O}', '\mu_{Si}', '\mu_{H}', '\mu_{Fe}', '\mu_{Sn}', '\mu_{Cr}', '\mu_{Ti}', '\mu_{Ir}', '\mu_{Mg}', '\mu_{Ca}', '\mu_{Zn}', '\mu_{Co}', '\mu_{Zr}', '\mu_{Hf}', '\mu_{Ta}', '\mu_{Ge}');
    xlabel('Equilibrium T (K)')
    ylabel('Chem Potentials (eV)')
    ylim([-12 0])




    %% plot the stoichiometry from equilibrium
    figure(7)
    [Ga_O_stoich, num_Ga, num_O, num_Si, num_H, num_Fe, num_Sn, num_Cr, num_Ti, num_Ir, num_Mg, num_Ca, num_Zn, num_Co, num_Zr, num_Hf, num_Ta, num_Ge] =  Ga2O3_stoich(equilib_dark_sol, conditions, defects);

    f7line(1) = plot(conditions.T_equilibrium,log10(num_Ga),'DisplayName','Ga');
    f7line(2) = plot(conditions.T_equilibrium,log10(num_O),'DisplayName','O');
    f7line(3) = plot(conditions.T_equilibrium,log10(num_Si),'DisplayName','Si');
    f7line(4) = plot(conditions.T_equilibrium,log10(num_H),'DisplayName','H');
    f7line(5) = plot(conditions.T_equilibrium,log10(num_Fe),'DisplayName','Fe');
    f7line(6) = plot(conditions.T_equilibrium,log10(num_Sn),'DisplayName','Sn');
    f7line(7) = plot(conditions.T_equilibrium,log10(num_Cr),'DisplayName','Cr');
    f7line(8) = plot(conditions.T_equilibrium,log10(num_Ti),'DisplayName','Ti');
    f7line(9) = plot(conditions.T_equilibrium,log10(num_Ir),'DisplayName','Ir');
    f7line(10) = plot(conditions.T_equilibrium,log10(num_Mg),'DisplayName','Mg');
    f7line(11) = plot(conditions.T_equilibrium,log10(num_Ca),'DisplayName','Ca');
    f7line(12) = plot(conditions.T_equilibrium,log10(num_Zn),'DisplayName','Zn');
    f7line(13) = plot(conditions.T_equilibrium,log10(num_Co),'DisplayName','Co');
    f7line(14) = plot(conditions.T_equilibrium,log10(num_Zr),'DisplayName','Zr');
    f7line(15) = plot(conditions.T_equilibrium,log10(num_Hf),'DisplayName','Hf');
    f7line(16) = plot(conditions.T_equilibrium,log10(num_Ta),'DisplayName','Ta');
    f7line(17) = plot(conditions.T_equilibrium,log10(num_Ge),'DisplayName','Ge');

    for i = 1:17
        f7line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f7line(i).DisplayName},size(f7line(i).XData)));
    end

    title('Atomic Concentrations')
    legend('Ga','O','Si','H','Fe','Sn','Cr', 'Ti', 'Ir', 'Mg', 'Ca', 'Zn', 'Co', 'Zr', 'Hf', 'Ta', 'Ge')
    xlabel('Equilibrium T (K)')
    ylabel('Log_{10} Concentrations (#/cm3)')
    ylim([conditions.plot_log10_min conditions.plot_log10_max])





    % analyze the distribution of each element across defects

    % here just using it for the frozen ones but could extend this for all
    % the elements inclduing those with fixed mu

    for ii = 1:numel(conditions.indices_of_fixed_elements)

        figure(7+ii)
        clf
        hold on
        [defects_with_element_fraction,~, defects_with_element_names, ~] = contains_element(equilib_dark_sol, defects, conditions.indices_of_fixed_elements(ii));

        for i=1:size(defects_with_element_fraction,2)
            plot(equilib_dark_sol.T_equilibrium,defects_with_element_fraction(:,i))
            %     eval(strcat("f",int2str(7+ii),"line(i) = plot(equilib_dark_sol.T_equilibrium,defects_with_element_fraction(:,i),'DisplayName',defects_with_element_names(i));"))
            %     eval(strcat("f",int2str(7+ii),"line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f",int2str(7+ii),"line(i).DisplayName},size(f"int2str(7+ii),"line(i).XData)));"));
        end

        title('Fraction of Element in Defects')
        % legend(defects_with_element_names')
        xlabel('Equilibrium T (K)')
        ylabel('Fraction of Element in Each Defect')
        ylim([0 1])

    end









    %%% clean up memory
    clear datacursormode f1line f11line f2line f3line f6line f7line i
    %%%%%%%%%%%%% end Tasks 1 and 2 (equilib and fullquenching)












% end  % end the for loop with index "index"
































% % %%%%%%%%%%%%%%% Task Trange 3 - quenching but setting interstitials generated at the equilibrium temp to zero   %%%%%%%%%%%%%%%
% %
% % % take out the interstitials by hand.  These defect and chargestate numbers
% % % only work for the full Varley native defect set (not validation)
% %
% % equilibrium_minus_interstitials_dark_sol = equilib_dark_sol;
% %
% % if strcmp(defect_set_flag, 'validation')
% %     %% If using the validation data set the interstitials are column 2 and 4.
% %     equilibrium_minus_interstitials_dark_sol.defects(:,2)=0;    % Ga_i
% %     equilibrium_minus_interstitials_dark_sol.chargestates(:,6:8)=0;
% %     equilibrium_minus_interstitials_dark_sol.defects(:,4)=0;   % O_i
% %     equilibrium_minus_interstitials_dark_sol.chargestates(:,12:14)=0;
% % elseif strcmp(defect_set_flag, 'full_Varley_set')
% %     %% if using the full Varley set, the interstitials are defects numbers 6, 13, and 14, and charge states 26,27,28 (Gai), and 52-56 (Oi and Osi)
% %     equilibrium_minus_interstitials_dark_sol.defects(:,6)=0;    % Ga_i
% %     equilibrium_minus_interstitials_dark_sol.chargestates(:,26:28)=0;
% %     equilibrium_minus_interstitials_dark_sol.defects(:,13)=0;   % O_i
% %     equilibrium_minus_interstitials_dark_sol.chargestates(:,52:54)=0;
% %     equilibrium_minus_interstitials_dark_sol.defects(:,14)=0;   % O_si
% %     equilibrium_minus_interstitials_dark_sol.chargestates(:,55:56)=0;
% % end
% %
% % tic
% % [quench_minus_interstitials_dark_sol] = defect_fullquench_dark(equilibrium_minus_interstitials_dark_sol, conditions, defects);
% % toc
% %
% % quench_minus_interstitials_dark_sol_out = cat(2,quench_minus_interstitials_dark_sol.T_equilibrium, quench_minus_interstitials_dark_sol.Nd, quench_minus_interstitials_dark_sol.Na, quench_minus_interstitials_dark_sol.EFn, quench_minus_interstitials_dark_sol.EFp, quench_minus_interstitials_dark_sol.n, quench_minus_interstitials_dark_sol.p, quench_minus_interstitials_dark_sol.defects);
% % save_name = strcat(conditions.name_root,'_quench_minus_interstitials_sol.dat');
% % writematrix(quench_minus_interstitials_dark_sol_out,save_name,'Delimiter','tab','WriteMode','overwrite')
% % disp(strcat(save_name,' written.'))
% %
% %
% % %% plot total defect concentrations
% % figure(6)
% % f6line(1) = plot(conditions.T_equilibrium,log10(quench_minus_interstitials_dark_sol.n),'ro');
% % f6line(2) = plot(conditions.T_equilibrium,log10(quench_minus_interstitials_dark_sol.p),'bo');
% % f6line(3) = plot(conditions.T_equilibrium,log10(abs(quench_minus_interstitials_dark_sol.Nd - quench_minus_interstitials_dark_sol.Na)),'k-');
% % for i = 1:3
% %     f6line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f6line(i).DisplayName},size(f6line(i).XData)));
% % end
% % legend('n','p','|Nd-Na|')
% % for i=1:defects.num_defects
% %     f6line(i+3) = plot(conditions.T_equilibrium,log10(quench_minus_interstitials_dark_sol.defects(:,i)),'DisplayName', defects.defect_names(i));
% %     f6line(i+3).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f6line(i+3).DisplayName},size(f6line(i+3).XData)));
% % end
% % legend show
% % title('Concentrations for Quenching all but O_i and Ga_i to 300 K')
% % xlabel('T for Equilibriation (K)')
% % ylabel('Log_{10} Concentrations at 300 K (#/cm3)')
% % ylim([10 23])
% % datacursormode.Enable='on';   %turn on the datatip to give name of defect on mouse hover
% %
% %
% %
% % %% add quenching minus interstitials results to Fig. 4
% % figure(4)
% % plot(conditions.T_equilibrium,quench_minus_interstitials_dark_sol.EFp,'k-')
% % plot(conditions.T_equilibrium,quench_minus_interstitials_dark_sol.EFn,'g-')
% % legend('Ec(T)','Ec(T_{quench})','Ev(T)','EFp(T) Equilib','EFn(T) Equilib', 'EFp(T) Quench', 'EFn(T) Quench', 'EFp(T) Quench, no interstitials', 'EFn(T) Quench, no interstitials')
% % xlabel('Equilibrium T (K)')
% % ylabel('Ec, Ev, EFn, EFp (eV)')
%
%
%
%
% %%%%% the tasks below here reqire a single temperature
% if size(conditions.T_equilibrium,2)~=1
%     disp('calculation of defects vs doping at fixed T requires only one T be in conditions.T_equilibrium - chose one T rather than a vector of T values')
%     yn = input('do you want to run the other tasks at a fixed temperature (must enter string y or n)');
%
%     if isempty(yn)
%         error('stopping... other parameter sweep tasks will be skipped')
%
%     elseif strcmp(yn,'y')
%         conditions.T_equilibrium = input('Enter one fixed temperature (K) for the rest of the calculations.  default T is set to 1300 K if you hit enter');
%         if isempty(conditions.T_equilibrium)
%             conditions.T_equilibrium = 1300;
%             for i=1:10
%                 disp('T set to 1000 C or 1300 K automatically')
%             end
%         else
%             disp('calculation of defects vs doping at fixed T requires only one T be in conditions.T_equilibrium - chose one T rather than a vector of T values')
%         end
%     else
%         error('stopping... other parameter sweep tasks will be skipped')
%     end
% end
%
%
%
%
% %%% Task Doping range 1:  explore defects vs doping at constant T %%%%%
%
%
% % make a temp holder for the shallow doping set at the start of the script,
% % reset back to that value once done.
% temp_Nd = conditions.Nd;
% temp_Na = conditions.Na;
%
% log_doping = 10:0.1:20;  %log10 of the doping range you want to sweep over
% num_dope = size(log_doping,2);
% defect_holder = zeros(num_dope,defects.num_defects);
% % chargestate_holder = zeros(num_dope,defects.num_chargestates);
% atom_holder = zeros(num_dope,6);    % change the 5 to a 6 once Sn added
% equilib_carrier_holder = zeros(num_dope,3);
% quench_carrier_holder = zeros(num_dope,3);
%
% for i = 1:num_dope
%     conditions.Nd = 10^log_doping(i);  %vary n-type shallow doping
% %     conditions.Na = 10^log_doping(i);  %vary p-type shallow doping
%     [equilib_dark_sol] = defect_equilibrium_dark_with_frozen(conditions, defects);
%     defect_holder(i,:) = equilib_dark_sol.defects;
% %     chargestate_holder(i,:) = equilib_dark_sol.defects;
%     equilib_carrier_holder(i,1) = equilib_dark_sol.n;
%     equilib_carrier_holder(i,2) = equilib_dark_sol.p;
%     equilib_carrier_holder(i,3) = abs(equilib_dark_sol.Nd - equilib_dark_sol.Na);  %no need to duplicate this for quenching
%     [fullquench_dark_sol] = defect_fullquench_dark(equilib_dark_sol, conditions, defects);
%     quench_carrier_holder(i,1) = fullquench_dark_sol.n;
%     quench_carrier_holder(i,2) = fullquench_dark_sol.p;
%     [Ga_O_stoich, atom_holder(i,1), atom_holder(i,2), atom_holder(i,3), atom_holder(i,4), atom_holder(i,5), atom_holder(i,6)] =  Ga2O3_stoich(equilib_dark_sol, conditions, defects);
% end
%
%
% figure(7)
% f7line(1) = plot(log_doping,log10(equilib_carrier_holder(:,1)),'ro','DisplayName','n_equilib');
% f7line(2) = plot(log_doping,log10(equilib_carrier_holder(:,2)),'bo','DisplayName','p_equilib');
% f7line(3) = plot(log_doping,log10(equilib_carrier_holder(:,3)),'k-','DisplayName','|Nd-Na|');  %this doesnt change for quenching so dont need to replot
% f7line(4) = plot(log_doping,log10(quench_carrier_holder(:,1)),'rx','DisplayName','n_quech');
% f7line(5) = plot(log_doping,log10(quench_carrier_holder(:,2)),'bx','DisplayName','p_quench');
% % for i = 1:5
% %     f7line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f7line(i).DisplayName},size(f7line(i).XData)));
% % end
% legend('n_{equilib}','p_{equilib}','|Nd-Na|','n_{quench}','p_{quench}')
%
% for i=1:defects.num_defects  %defects dont change for quenching so only need to plot once
%     f7line(i+5) = plot(log_doping,log10(defect_holder(:,i)),'DisplayName', defects.defect_names(i));
%     f7line(i+5).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f7line(i+5).DisplayName},size(f7line(i+5).XData)));   %the i+5 here is becasue we have 5 lines above for carrier densities
% end
% legend show
% title('Concentrations at Fixed T vs Shallow Doping')
% xlabel('Log_{10}(Nd) (#/cm^3)')
% ylabel('Log_{10} Concentrations (#/cm^3)')
% ylim([10 23])
% datacursormode.Enable='on';   %turn on the datatip to give name of defect on mouse hover
%
%
%
% figure(8)
% f8line(1) = plot(log_doping,log10(equilib_carrier_holder(:,1)),'ro','DisplayName','n_equilib');
% f8line(2) = plot(log_doping,log10(equilib_carrier_holder(:,2)),'bo','DisplayName','p_equilib');
% f8line(3) = plot(log_doping,log10(equilib_carrier_holder(:,3)),'k-','DisplayName','|Nd-Na|');  % this doesnt change for quenching so dont need to replot
% f8line(4) = plot(log_doping,log10(quench_carrier_holder(:,1)),'rx','DisplayName','n_quech');
% f8line(5) = plot(log_doping,log10(quench_carrier_holder(:,2)),'bx','DisplayName','p_quench');
% for i = 1:5
%     f8line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f8line(i).DisplayName},size(f8line(i).XData)));
% end
% legend('n_{equilib}','p_{equilib}','|Nd-Na|','n_{quench}','p_{quench}')
%
% for i=1:6  % defects dont change for quenching so only need to plot once.  Plot the atomic concentrations
%     f8line(i+5) = plot(log_doping,log10(atom_holder(:,i)),'DisplayName', defects.elementnames(i));
%     f8line(i+5).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f8line(i+5).DisplayName},size(f8line(i+5).XData)));   %the i+5 here is becasue we have 5 lines above for carrier densities
% end
% legend show
% title('Concentrations at Fixed T vs Shallow Doping')
% xlabel('log_{10}(Nd) (#/cm^3)')
% ylabel('Log_{10} Concentrations (#/cm^3)')
% ylim([10 23])
% datacursormode.Enable='on';   %turn on the datatip to give name of defect on mouse hover
%
% %% reset the Nd and Na now that we're done
% conditions.Nd = temp_Nd;
% conditions.Na = temp_Na;
% clear temp_Na temp_Nd f8line f7line datacursormode
% %%%%%%%%%%%%%%% end of Task 4 loop over doping
% %%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%
% % % %  Task 5 loop over mu for some element
% %% to use this, you must select the right line in the for loop below for the selected element
%
% for i=1:10
%     disp('dont forget to reset the conditions.muT_equilibrium values after running task 5')
% end
%
% %%%%%%%%% take the mu set at top of script, reset it back after this task
% temp_mu = conditions.muT_equilibrium;
%
% if size(conditions.T_equilibrium,2)~=1
%     error('Task 5 calculation of defects vs mu at fixed T requires only one T be in conditions.T_equilibrium - chose one T rather than a vector of T values')
% end
%
% mu_vec = -6:0.1:0;  % enter the range to search over
% num_mu_in_mu_vec = size(mu_vec,2);
% % chargestate_holder = zeros(num_mu_in_mu_vec,defects.num_chargestates);
% defect_holder = zeros(num_mu_in_mu_vec,defects.num_defects);
% atom_holder = zeros(num_mu_in_mu_vec,defects.numelements);    % 6 corresponds to the list Ga O Si H Fe Sn
% equilib_carrier_holder = zeros(num_mu_in_mu_vec,3);
% quench_carrier_holder = zeros(num_mu_in_mu_vec,3);
%
% for i = 1:num_mu_in_mu_vec
%     conditions.muT_equilibrium(:,conditions.Task5_element_num_to_vary_mu) = mu_vec(i);   % assign the desired value of mu to the right element
%     [equilib_dark_sol] = defect_equilibrium_dark_with_frozen(conditions, defects);  % run the calc
%     defect_holder(i,:) = equilib_dark_sol.defects;
% %     chargestate_holder(i,:) = equilib_dark_sol.defects;
%     equilib_carrier_holder(i,1) = equilib_dark_sol.n;
%     equilib_carrier_holder(i,2) = equilib_dark_sol.p;
%     equilib_carrier_holder(i,3) = abs(equilib_dark_sol.Nd - equilib_dark_sol.Na);
%     [fullquench_dark_sol] = defect_fullquench_dark(equilib_dark_sol, conditions, defects);
%     quench_carrier_holder(i,1) = fullquench_dark_sol.n;
%     quench_carrier_holder(i,2) = fullquench_dark_sol.p;
%     quench_carrier_holder(i,3) = abs(fullquench_dark_sol.Nd - fullquench_dark_sol.Na);
%     [Ga_O_stoich, atom_holder(i,1), atom_holder(i,2), atom_holder(i,3), atom_holder(i,4), atom_holder(i,5), atom_holder(i,6)] =  Ga2O3_stoich(equilib_dark_sol, conditions, defects);
% end
%
% figure(9)
% f9line(1) = plot(mu_vec,log10(equilib_carrier_holder(:,1)),'ro','DisplayName','n');
% f9line(2) = plot(mu_vec,log10(equilib_carrier_holder(:,2)),'bo','DisplayName','p');
% f9line(3) = plot(mu_vec,log10(equilib_carrier_holder(:,3)),'k-','DisplayName','|Nd-Na|');
% f9line(4) = plot(mu_vec,log10(quench_carrier_holder(:,1)),'rx','DisplayName','n_quech');
% f9line(5) = plot(mu_vec,log10(quench_carrier_holder(:,2)),'bx','DisplayName','p_quench');
% for i = 1:5
%     f9line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f9line(i).DisplayName},size(f9line(i).XData)));
% end
% legend('n_{equilib}','p_{equilib}','|Nd-Na|','n_{quench}','p_{quench}')
% for i=1:defects.num_defects
%     f9line(i+5) = plot(mu_vec,log10(defect_holder(:,i)),'DisplayName', defects.defect_names(i));
%     f9line(i+5).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f9line(i+3).DisplayName},size(f9line(i+3).XData)));
% end
% legend show
% title(strcat('Concentrations at Fixed T vs_',defects.elementnames(conditions.Task5_element_num_to_vary_mu),' Chemical Potential'))
% xlabel(strcat(defects.elementnames(conditions.Task5_element_num_to_vary_mu),' Chemical potential (eV)'))
% ylabel('Log_{10} Concentrations (#/cm^3)')
% ylim([10 23])
% datacursormode.Enable='on';   %turn on the datatip to give name of defect on mouse hover
%
%
%
%
% figure(10)
% f10line(1) = plot(mu_vec,log10(equilib_carrier_holder(:,1)),'ro','DisplayName','n_equilib');
% f10line(2) = plot(mu_vec,log10(equilib_carrier_holder(:,2)),'bo','DisplayName','p_equilib');
% f10line(3) = plot(mu_vec,log10(equilib_carrier_holder(:,3)),'k-','DisplayName','|Nd-Na|');  % this doesnt change for quenching so dont need to replot
% f10line(4) = plot(mu_vec,log10(quench_carrier_holder(:,1)),'rx','DisplayName','n_quech');
% f10line(5) = plot(mu_vec,log10(quench_carrier_holder(:,2)),'bx','DisplayName','p_quench');
% for i = 1:5
%     f10line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f10line(i).DisplayName},size(f10line(i).XData)));
% end
% legend('n_{equilib}','p_{equilib}','|Nd-Na|','n_{quench}','p_{quench}')
%
% for i=1:6  % defects dont change for quenching so only need to plot once.  Plot the atomic concentrations
%     f10line(i+5) = plot(mu_vec,log10(atom_holder(:,i)),'DisplayName', defects.elementnames(i));
%     f10line(i+5).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f10line(i+5).DisplayName},size(f10line(i+5).XData)));   %the i+5 here is becasue we have 5 lines above for carrier densities
% end
% legend show
% title('Element Concentrations at Fixed T vs Shallow Doping')
% xlabel(strcat(defects.elementnames(conditions.Task5_element_num_to_vary_mu),' Chemical potential (eV)'))
% ylabel('Log_{10} Concentrations (#/cm^3)')
% ylim([10 23])
% datacursormode.Enable='on';   %turn on the datatip to give name of defect on mouse hover
%
% conditions.muT_equilibrium = temp_mu ;   % reset the mu values back to what they were before running task 5
% clear temp_mu f9line f10line datacursormode
%
% %%%%%%%%%%%%%%%%%% end of Task 5 %%%%%%%%%%%%%%



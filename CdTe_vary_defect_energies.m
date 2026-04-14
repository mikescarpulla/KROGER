%% clear workspace and figures   %%%%%%%%%%%%
clear all
clc

%% set Physical Constants  %%%%%%%%%
conditions.q = 1.602176565e-19;
conditions.h = 6.62606957e-34;
conditions.kB = 8.6173324e-5;
conditions.mo = 9.1093837e-31;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% Cd and Te compositions for the boundaries of the solidus region and the mu values for liquidus
% hard code the excess Cd or Te points here
Cd_rich_T = [700 840 900 1000 1100]; % list of measured points
Cd_rich_solidus = [1e17 1e18 1e19 1e18 1e17]; %
Cd_rich_liquidus_mu = [0 -1; -1.3 -1.4; 0 -1; -1.3 -1.4; -1.5 1];  %muCd = col1 muTe=col2

Te_rich_T = [750 845 900 950 1045];
Te_rich_solidus = [1e17 1e18 1e19 1e20 1e19];
Te_rich_liquidus_mu = [0 -1; -1.3 -1.4; 0 -1; -1.3 -1.4; -1.5 1];  %muCd = col1 muTe=col2

excess_ratio = 5; % cutoff ratios to keep scenarios

conditions.save_fname = "CdTe_excess_sweep_dHo.txt";


%%  Select the .mat file database of defect properties and load it %%%%%%%%
%%% note: defects and chargestates listed down, one row per.  Elements
%%% listed as columns in variables like temp_defects.cs_dm.  These conventions must
%%% be obeyed so that the matrix and vector multiplicaitons used in calcs
%%% work out right - e.g. when mu vectors and dm vectors are multiplied to
%%% get formation energies for each defect.

conditions.defect_db_fname = [];
conditions.defect_db_pname = [];
[conditions.defect_db_fname,conditions.defect_db_pname] = uigetfile('*.mat','Select the defect database file you want to use');

% if ~ischar(conditions.defect_db_fname) || isempty(conditions.defect_db_fname)   %if the user doesnt pick one in the UIgetfile, do it manually here
%     conditions.defect_db_fname = 'CdTe_defects_Intuon_10062025.mat';
% end

if ~ischar(conditions.defect_db_pname) || isempty(conditions.defect_db_pname)
    conditions.defect_db_pname = strcat(pwd,'\');  % the pwd command leaves off the last \ so have to add it.
end

load(strcat(conditions.defect_db_pname,conditions.defect_db_fname))
% defect_db_name_root = conditions.defect_db_fname(1:(size(conditions.defect_db_fname,2) - 4)    );  % this takes the filename and strips the . and extension - keep it until savenames established

if numel(fieldnames(defects)) < 18
    error('Defects database file seems to be missing some required parts - carefully check it and try again')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conditions.save_pname = conditions.defect_db_pname;
clear save_fname_root save_fname_root_length defect_set_flag defect_db_name_root
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% Calc and Display Options  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

conditions.plot_log10_min = 12;   %%% set a threshold - if a quantity doesn't go above this at some temperature then dont bother plotting it
conditions.plot_log10_max = 20;   %%% set upper bound for plots too.  The real numbers are 10^x for these x's

% save to file options
% min cutoff concentration for significant concentrations
conditions.save_min = 1e12;  % set a minimum threshold for a defect or chargestate to be saved to the reduced output files (real number not log of it)

% define the material
conditions.stoich_flag = 'CdTe';

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
        disp(msg)
    end
elseif conditions.num_fixed_defects==0
    disp('No Defects fixed.  Proceeding')
end

%%%  end section on fixed defects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% save a copy of defects and of conditions
temp_defects = defects;
temp_conditions = conditions;



counter = 0;
for i=1:2  % loop over Cd rich or Te rich case with i

    if i==1
        temp_conditions.T_equilibrium = Cd_rich_T;
        temp_conditions.muT_equilibrium = Cd_rich_liquidus_mu;
        temp_conditions.num_T_equilibrium = numel(Cd_rich_T);
    elseif i==2
        temp_conditions.T_equilibrium = Te_rich_T;
        temp_conditions.muT_equilibrium = Te_rich_liquidus_mu;
        temp_conditions.num_T_equilibrium = numel(Te_rich_T);
    end



%% now we have the T vector we are going to use for this scenario



    %%%%% Determine how to set pressures, chemical potentials, concentrations.. of the elements

    %% set total pressure and units %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    temp_conditions.P_tot = 1;   % total pressure in atm
    temp_conditions.P_units = 'atm';

    %% Set up the elements to use in the model %%%%%%%%%%%%%
    if ~isfield(conditions,'num_elements') || sum(isempty(temp_conditions.num_elements))
        temp_conditions.num_elements = size(defects.cs_dm,2);
    end

    %% Set the scenario used to determine the chemical potential of the host material.
    %%% different options for setting the host material's chemical potentials according to equilibrium equations or other

    % for T-independent mu conditions (set one mu like mu=0 for rich/poor conditions)
    temp_conditions.mu_CdTe = -1.39;        % formation energy of CdTe from Intuon

    %% set chemical potentials depending on the method for each element
    temp_conditions.mu_constant_flag = ones(1,temp_conditions.num_elements);
    temp_conditions.mu_from_2nd_phases_flag = zeros(1,temp_conditions.num_elements);   % set all to zero/false first then overwrite if needed
    temp_conditions.muT_equilibrium = -30*ones(temp_conditions.num_T_equilibrium, defects.num_elements); %


    % logic checks
    mu_set_exclusively_by_one_method_check = ((temp_conditions.mu_constant_flag + temp_conditions.mu_from_2nd_phases_flag)-1) ~=0 ;
    temp_conditions.mu_set_flag = (temp_conditions.mu_constant_flag + temp_conditions.mu_from_2nd_phases_flag)>0;

    if any(mu_set_exclusively_by_one_method_check)
        disp(mu_set_exclusively_by_one_method_check)
        error('mu for each element has to be set either as a constant or by a phase boundary.  The line above has a 1 in the position of the element causing the problem')
    end

    temp_conditions.fixed_conc_flag = zeros(1,temp_conditions.num_elements);
    temp_conditions.fixed_conc_values = zeros(1,temp_conditions.num_elements);  % this will hold target values of elements
    temp_conditions.fixed_elements_mu_ranges = zeros(temp_conditions.num_elements,2);   % create  matrix of right size for ALL the elements - only will pay attention to those which are fixed

    % default is no elements fixed - catch if nothing entred.
    if ~isfield(conditions,'fixed_conc_flag') || sum(isempty(temp_conditions.fixed_conc_flag))
        temp_conditions.fixed_conc_flag = zeros(defects.num_elements,1);
    elseif ~isfield(conditions,'fixed_conc_values') || sum(isempty(temp_conditions.fixed_conc_values))
        temp_conditions.fixed_conc_values = zeros(defects.num_elements,1);
    end
    temp_conditions.fixed_conc_values = temp_conditions.fixed_conc_values.*temp_conditions.fixed_conc_flag;
    temp_conditions.indices_of_fixed_elements = find(temp_conditions.fixed_conc_flag);
    temp_conditions.num_fixed_elements = numel(temp_conditions.indices_of_fixed_elements);
    if temp_conditions.num_fixed_elements~=0
        for z=1:temp_conditions.num_fixed_elements
            msg = strcat('WARNING: Fixed concentration of element used for...',defects.element_names(temp_conditions.indices_of_fixed_elements(z)),'...Be sure you meant to do this.');
            disp(msg)
        end
    elseif temp_conditions.num_fixed_elements==0
        disp('No Eelements have fixed concentrations.  Proceeding')
    end
    clear msg impurity_mu_boundary_limits mu_set_exclusively_by_one_method_check max_mu_range start_kBT fine_factor
    clear mu_Cd mu_Te
    %%% end section on fixed elements %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    %% logic checks on mu or concentraitons set for each element
    % catch it if by accident the mu and defect conc are both set for any
    % elements

    neither_mu_nor_conc_flag = (temp_conditions.mu_set_flag .* temp_conditions.fixed_conc_flag);
    both_mu_and_conc_flag = temp_conditions.mu_set_flag + temp_conditions.fixed_conc_flag ==0;

    if any(neither_mu_nor_conc_flag)
        disp(neither_mu_nor_conc_flag)
        error('Must set either the mu or concentation for each element.')
    elseif any(both_mu_and_conc_flag)
        disp(both_mu_and_conc_flag)
        error('Must set either the mu or the concentation for each element, not both.')
    elseif size(temp_conditions.muT_equilibrium,2)~=temp_conditions.num_elements
        error('number of columns in mu array not equal to number of elements declared')
    elseif sum(sum(isinf(temp_conditions.muT_equilibrium)))
        error('Infinite mu detected - this happens when there is no model for G0 for the T requested for some substance(s)')
    end

    clear both_mu_and_conc_flag neither_mu_nor_conc_flag
    %%% end handling matrix chemical potentials %%%%%%%%%%%%


    % this has to come after the temperature vector is set
    %%% Ga2O3 %%

    temp_conditions.vibent_T0 = 200; % This is the characteristic To for the average phonon mode in Ga2O3, as determined from the wo that makes the Debeye temperature come out to about 740 K based on the Cv(T) data found online.  To where x = hbar*w0/kBT = T0/T, wherein w0 is the mean phonon energy from the Debeye temperature determined from Cv(T).

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

    % these take the equilibrium T values from conditions and generate the
    % needed T-dependent values.  The fraction of Eg change caused by VB vs CB
    % can be modified above.
    temp_conditions.NcT_equilibrium = temp_conditions.NcRef*(temp_conditions.T_equilibrium/temp_conditions.TRef).^1.5;
    temp_conditions.NvT_equilibrium = temp_conditions.NvRef*(temp_conditions.T_equilibrium/temp_conditions.TRef).^1.5;
    delta_EgT_equilibrium = (temp_conditions.varshini_a .* temp_conditions.T_equilibrium.^2)./(temp_conditions.T_equilibrium + temp_conditions.varshini_b);  %delta is a positive number
    temp_conditions.EgT_equilibrium = temp_conditions.Eg0 - delta_EgT_equilibrium;
    temp_conditions.EcT_equilibrium = temp_conditions.Eg0*ones(size(temp_conditions.EgT_equilibrium)) - temp_conditions.EcT_fraction*delta_EgT_equilibrium;
    temp_conditions.EvT_equilibrium = zeros(size(temp_conditions.EgT_equilibrium)) + temp_conditions.EvT_fraction*delta_EgT_equilibrium;
    clear delta_EgT_equilibrium delta_EgT_fullquench   %clean up memory after using these

    % site densities in the Ga2O3 lattice (num FU/unit cell / vol for unit cell)
    N_Cd = 1.47e22;    %site 1
    N_Te = 1.47e22;  % site 2
    N_iCd = 1.47e22;  %site 3
    N_iTe = 1.47e22;   %site 4
    temp_conditions.num_sites = [N_Cd; N_Te; N_iCd; N_iTe];   % note this needs to be a column vector, dont change it to a row vector

    % % Calculate the numerical prefactor for each defect from defects.degen_factor_config, defects.degen_factor_elec, defects.cs_num_each_site, and temp_conditions.num_sites

    %% defects.cs_num_each_site tells you the number of primitive unit cells needed to form each defect (assuming the sites of the crystal, including distinct interstitials like the ia, ib, ... in b-ga2O3, are counted such that there is one per primitive unit cell
    defects.cs_site_prefactor = sort((ones(defects.num_cs,1)*temp_conditions.num_sites')./defects.cs_num_each_site,2);      % calculates the max number of each chargestate that could be formed, given the availability of each site type in the crystal (or unit cells, or supercells, whatever the basis for counting).  Then sorts by the rows ascending, so the first column becomes the max number we could form.  Its done with a sort becasue you end up with Inf when you divide by sites you need zero of
    defects.cs_site_prefactor = defects.cs_site_prefactor(:,1);  % toss out the other columns keeping only the first one, which should be the limiting one
    defects.cs_tot_prefactorfactor = defects.cs_degen_factor_config .* defects.cs_degen_factor_elec .* defects.cs_site_prefactor;

    clear N_Cd N_Te N_iCd N_iTe
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








    % vector of energies to add to the numbers in the defect database (can be negative or positive)
    % these will be permuted across the different defects
    d_nrg = -1 : 0.5 : 2;
    num_nrgs = numel(d_nrg);

    % total number of scenarios we want to check for each T
    num_scenarios = num_nrgs^temp_defects.num_defects;

    if num_scenarios > 5e6
        error("number of scenarios >5e6 - matrices will be too big")
    end

    data_holder = zeros(num_scenarios, 2 + defects.num_cs);

    % we will do this by using ndgrid to construct meshed arrays, then index
    % into them using linear indexing.  what we will do is build a variabe command string.  Rather than making a
    % fixed nested for loop, we can build the right number of permutations up front and then use the trick of linear indexing into multidimensional arrays to run the whole thing as one for loop.
    % So build up a list of permuted scenarios up front, then just go through and run each one.


    % numbers of sites per cm3 in perfect lattice
    Nsites_Cd = temp_conditions.num_sites(1);
    Nsites_Te = temp_conditions.num_sites(2);
    Num_Cd = Nsites_Cd;
    Num_Te = Nsites_Te;


    max_lim = excess_ratio * max([Cd_rich_solidus Te_rich_solidus],[],"all"); % only accept solutions where the excess atom numbers are within a ratio of the min and max measured points
    min_lim = 1/excess_ratio * min([Cd_rich_solidus Te_rich_solidus],[],"all"); % same for min




    if i==1
        temp_conditions.T_equilibrium = Cd_rich_T;
        excess_vec = Cd_rich_solidus;
    elseif i==2
        temp_conditions.T_equilibrium = Te_rich_T;
        excess_vec = Cd_rich_solidus;
    end

    temp_conditions.num_T_equilibrium = numel(temp_conditions.T_equilibrium);
    temp_conditions.kBT_equilibrium = temp_conditions.kB * temp_conditions.T_equilibrium;

    str1 = "";
    str2 = "";
    for j = 1:temp_defects.num_defects
        str1 = str1 + "deltaE_" + int2str(j) + " ";
        str2= str2 + ", d_nrg";
    end
    str2 = char(str2); % strip off the first comma
    str2 = string(str2(3:end));
    text_command = "[" + str1 + "] = ndgrid(" + str2 + ");" ;
    eval(text_command)

    % we have now one variable per defect called deltaE_x where x is the
    % defect number.


    % at this point, we have generated a hyper-array (dimensions = num_defects) whcih holds permuted combinations of all the values in the d_nrg vector
    % across all of the defects.  we can use linear indexing to get the right values for all of the defects




    for l = 1:num_scenarios  % we can index in to a 3x3 array 1 = (1,1,1), 2=(1,1,2),.. 4 = (1,2,1), 5=(1,2,2) etc - we jsut do it in N dimensions now.  It doesnt matter the order we check them as longas we get the corresponding dHo from all of the defects systematically
        for m = 1:temp_defects.num_defects
            for n = temp_defects.cs_indices_lo(m) : temp_defects.cs_indices_hi(m)   % index over the chargestates of the jth defect
                temp_defects.cs_Eform(n) = temp_defects.cs_Eform(n) + eval("deltaE_" + int2str(m) + "(" + int2str(l) + ")");  %find the dHo of the kth chargestate and add the jth delta_energ to it. Store the result in a temprary defects varaible that changes for each scenario
            end
        end
        % at this point we havechanged all the Eform values for this scenario

        % run the calculaiton at all the temperatures for Cd or Te rich
        % (i=1 or i=2 cases)
        [equilib_dark_sol] = defect_equilibrium_dark(temp_conditions, temp_defects);

        % calc stoichiometry
        [~,element_totals] = CdTe_stoich(equilib_dark_sol, temp_conditions, temp_defects);

        if i==1  % just pick out the Cd or Te depending on i and make it one column
            element_excess(:) = element_totals(:,1) - Nsites_Cd;
        elseif i==2
            element_excess(:) = element_totals(:,2) - Nsites_Te;
        end

        % check if this entire scenario is within the bounds (for all
        % temperatures)
        scenario_within_limits_flag = all(element_excess <= max_lim) && all(element_excess >= min_lim);

        if scenario_within_limits_flag
            counter = counter + 1
            chi2 = 1/temp_conditions.num_T_equilibrium * sum((excess_vec - element_excess).^2);
            data_line = [i chi2 temp_defects.cs_Eform'];
            data_holder(counter,:) = data_line;
        end
    end  %end l loop



end %end i loop


% write output file
headers = ['1=Cd_rich, 2=Te_rich' 'scenario chi2' temp_defects.cs_names'];
output_cell = [cellstr(headers); num2cell(data_holder)];
writecell(output_cell,strcat(temp_conditions.save_pname,temp_conditions.save_fname),'Delimiter','tab','WriteMode','overwrite');
disp('wrote element total concentrations to text file')


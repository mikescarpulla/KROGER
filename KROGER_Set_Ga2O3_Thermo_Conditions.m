
%% set Temperatures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vector of T values at which to calculate for equilibrium
% start at the highest T you want and then work down in order to fully
% exploit the adaptive guessing.  Also, for fixed elements the major time
% is spent finding the first guess, the optimization is fast by comparison once you fix 2 or more elements.
% So it's worth doing a fine T grid for long calcuations so you dont have
% to ever do it over.

conditions.T_equilibrium = 2100:-100:500;  % row vector of T in K at which to compute the equilirium

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
conditions.pO2 = 0.2;
% conditions.pO2 = 0.02;
% conditions.pO2 = 0.2*conditions.P_tot;  % for MOCVD
% conditions.pO2 = 1e-4;
% conditions.pO2 = 1e-6;
% conditions.pO2 = 1e-8;

% if using "pGa_vap-variable", need to also set a value
% conditions.pGa_vap = 1e-4;
% conditions.pGa_vap = 1e-8;
% conditions.pGa_vap = 1e-10;

% specify partial pressures of Ga2O or GaO
% conditions.pGa2O_vap = 1e-4;
% conditions.pGaO_vap = 1e-4;

if strcmp(conditions.mu_conditions_flag,'O2-1atm')   % equilibrium with O2 at pO2=ptot=1 atm.  This overrides anything you set accidentally above.
    conditions.pO2 = 1;
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
% 
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
% 1=Ga 2=O 3=Si 4=H 5=Fe 6=Sn 7=Cr 8=Ti 9=Ir 10=Mg 11=Ca 12=Zn 13=Co 14=Zr 15=Hf 16=Ta 17=Ge


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

conditions.fixed_conc_flag(6) = 1;    % fix Sn
conditions.mu_set_flag(6) = 0;              % flip the flag saying Sn is set by mu
conditions.fixed_conc_values(6) = 4.5e18;
% % 2e18 Sn is about -5.5 eV at 2100K  to -5 eV at 500 K
lo_mu_Sn = -6.5; % set range to search for the unknown mu
hi_mu_Sn = -3;
conditions.fixed_elements_mu_ranges(6,:) = [lo_mu_Sn hi_mu_Sn];

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
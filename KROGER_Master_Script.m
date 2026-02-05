% CdTe - these call scripts of the same name.  Edit them to control the desired calculation
KROGER_Load_CdTe_Database_and_Calc_Settings
KROGER_Set_CdTe_Thermo_Conditions  % this includes temperature vector so it needs to come before material and defects
KROGER_Set_CdTe_Material_Conditions
KROGER_Set_CdTe_Defect_Conditions

%% fix As concentration 
% % conditions.fixed_conc_values(5) = 1e17;     % specify the value
% % lo_mu_As = -10; % set range to search for the unknown mu
% % hi_mu_As = -2;
% % conditions.fixed_elements_mu_ranges(5,:) = [lo_mu_As hi_mu_As];

%% do the calc, post process the results, save and plot

% conditions.dopant_sweep = [1e15 1e16 1e17 1e18 1e19];  % As doping conc.
% conditions.dopant_sweep = 1e17;
% 
% conditions.num_dopant_sweep = numel(conditions.dopant_sweep);
% % 
% % for i = 1:conditions.num_dopant_sweep
% 
%     conditions.fixed_conc_values(5) = conditions.dopant_sweep(i); % element num 5 is As

    [equilib_dark_sol] = defect_equilibrium_dark(conditions, defects);
    % [equilib_dark_sol] = defect_equilibrium_dark_with_frozen(conditions, defects);

    KROGER_Save_Full_Equlibrium_Outputs

    [fullquench_dark_sol] = defect_fullquench_dark(equilib_dark_sol, conditions, defects);
    KROGER_Save_Full_Quench_Outputs

    % KROGER_grouping_defects

    % special_calc

    KROGER_Plot_Outputs
% end

% KROGER_Validate_and_Save_Conditions_and_Do_Calculations
% KROGER_Save_Full_Equlibrium_Outputs
% KROGER_Save_Full_Quench_Outputs
% KROGER_Plot_Outputs



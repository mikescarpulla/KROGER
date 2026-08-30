% CdTe - these call scripts of the same name.  Edit them to control the desired calculation
KROGER_Load_CdTe_Database_and_Calc_Settings
KROGER_Set_CdTe_Thermo_Conditions  % this includes temperature vector so it needs to come before material and defects
KROGER_Set_CdTe_Material_Conditions
KROGER_Set_CdTe_Defect_Conditions

% 
% % Ga2O3 - these call scripts of the same name.  Edit them to control the desired calculation
% KROGER_Load_Ga2O3_Database_and_Calc_Settings
% KROGER_Set_Ga2O3_Thermo_Conditions  % this includes temperature vector so it needs to come before material and defects
% KROGER_Set_Ga2O3_Material_Conditions
% KROGER_Set_Ga2O3_Defect_Conditions



% conditions.loop_over = 'Nd';
Nd_vec = [1e12 1e13 1e14 1e15 1e16 1e17 1e18 1e19 1e20];
num_Nd = numel(Nd_vec);

for script_i = 1:num_Nd

%% do the calc, post process the results, save and plot
KROGER_Validate_and_Save_Conditions_and_Do_Calculations

build up the solution variable here, where rows are for the looped variable

end

% KROGER_Save_Full_Equlibrium_Outputs
% KROGER_Save_Full_Quench_Outputs
% KROGER_Plot_Outputs
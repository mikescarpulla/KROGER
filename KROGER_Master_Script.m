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


%% do the calc, post process the results, save and plot
KROGER_Validate_and_Save_Conditions_and_Do_Calculations
KROGER_Save_Full_Equlibrium_Outputs
KROGER_Save_Full_Quench_Outputs
KROGER_Plot_Outputs

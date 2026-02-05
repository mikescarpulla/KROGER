% CdTe - these call scripts of the same name.  Edit them to control the desired calculation
KROGER_Load_CdTe_Database_and_Calc_Settings

% all probelms will start at some high T and go to a low T final state.
% The SQ calc will be done starting from the points specified here.  We use
% Tequlibrium to generate and hold all the needed things like mu, Eg, etc
% at 1 K spacings.  
conditions.SQ_T_start = 1500; % assign highest T to consider
conditions.SQ_T_end = 300; % end temp for cooling (usually 300 K)
conditions.SQ_T_step = 50;

KROGER_SQ_Set_CdTe_Thermo_Conditions  % this includes temperature vector so it needs to come before setting material and defects
% KROGER_SQ_Set_Ga2O3_Thermo_Conditions  % this includes temperature vector so it needs to come before setting material and defects

KROGER_Set_CdTe_Material_Conditions
% KROGER_Set_Ga2O3_Material_Conditions

KROGER_Set_CdTe_Defect_Conditions
% KROGER_Set_Ga2O3_Defect_Conditions




% now set things only found for SQ calc.  

% linear or exp cooling trajectory
conditions.SQ_lin_or_exp = 'lin';
% conditions.SQ_lin_or_exp = 'exp';


%% set the characteristic distance (equal distance in log scale)
chardist_decadal_spacing = 2;  % 1 means each decade, 0.1 means 1,2,3,4...(10 per decade)
min_chardist_10exponent = -4;
max_chardist_10exponent = 2;

%% Set cooling rate.  times in sec and rates in K/s
Trate_decadal_spacing = 2;  % 2 means step by 100's, 1 means each decade, 0.1 means 1,2,3,4...(10 per decade)
min_Trate_10exponent = -6;
max_Trate_10exponent = 0;


%% flesh things out
conditions.SQ_chardists = 10.^(min_chardist_10exponent : chardist_decadal_spacing : max_chardist_10exponent);
conditions.SQ_num_chardists = numel(conditions.SQ_chardists);
conditions.SQ_Trates = 10.^(min_Trate_10exponent : Trate_decadal_spacing : max_Trate_10exponent);
conditions.SQ_num_Trates = numel(conditions.SQ_Trates);
% now calculate cooltimes from Trates
if strcmp(conditions.SQ_lin_or_exp,'lin')
    conditions.SQ_cooltimes = (conditions.SQ_T_start-conditions.SQ_T_end)./conditions.SQ_Trates;  % if linear, this is the cooling time to go from Tmax-Tmin.
elseif strcmp(conditions.SQ_lin_or_exp,'exp')
    conditions.SQ_cooltimes = conditions.SQ_T_start./conditions.SQ_Trates;  % If exp, this is the characteristic decay time t0 in sec (time to decrease by 1/e) - can rewrite as exp(-t/to) with to=Tmax/rate
end
conditions.SQ_num_cooltimes = numel(conditions.SQ_Trates);




% actual function call 
[SQ_dark_sol, SQ_conditions] = SQ_calculation_diffusion(defects, conditions);

%% SQ modification by defect removal
modification.which_defects = [3 4];     % IDs of defects that need to modify
modification.new_value = [0 0];       

[mod_dark_sol] = modify_SQ_by_defect_removal(SQ_dark_sol, conditions, defects, modification);
%%

SQ_plot_outputs




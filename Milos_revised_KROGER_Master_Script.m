%% Instructions
% step 1 - fix everything in the set and load files
% step 2 - decide which calc to do (single, Na/Nd sweep, dopant sweep)
% step 3 - set saving terms if you want to auto save figures


% GaN - these call scripts of the same name.  Edit them to control the desired calculation
KROGER_Load_GaN_Database_and_Calc_Settings
KROGER_Set_GaN_Thermo_Conditions  % this includes temperature vector so it needs to come before material and defects
KROGER_Set_GaN_Material_Conditions
KROGER_Set_GaN_Defect_Conditions

%% Set controls (dopant sweep) Three calc commands (single, Na/Nd sweep, dopant sweep)

% masters.calc_command_flag = 'single';
% masters.calc_command_flag = 'Na/Nd sweep';
masters.calc_command_flag = 'dopant sweep';

% Single (Here you can set these often changed values instead of checking
% multiple scripts)
conditions.Nd = 0;
conditions.Na = 1e15;
conditions.mu_conditions_flag ='Ga-rich';

% Na/Nd sweep
masters.N_sweep = [1e15 1e16 1e17 1e18 1e19];
masters.num_N_sweep = numel(masters.N_sweep);
masters.N_sweep_flag = 'Na';
% masters.N_sweep_flag = 'Nd';

% dopant sweep (this does not substitute fixing the element in the
% thermo_conditions file)
masters.dopant_to_sweep_number = 4; % Ga N C Si Mg Ge O Be Zn Ca Cd H
masters.dopant_sweep = [1e15 1e16 1e17 1e18];
masters.num_dopant_sweep = numel(masters.dopant_sweep);

%% Set saving terms
masters.auto_save_flag = 'on';
masters.save_root = '/Users/mr.coolguymilo/Desktop/Scarpulla Lab Work/GaN Matlab Scripts/The ones Scarp revised/Calcs/New Master Script Testing/calc dopant sweep';
masters.auto_save = [2, 3];
%% Run Calcs
if (strcmp(masters.calc_command_flag, 'single'))
    calc_single(conditions, defects, masters)
elseif (strcmp(masters.calc_command_flag, 'Na/Nd sweep'))
    calc_N_sweep(conditions, defects, masters)
elseif (strcmp(masters.calc_command_flag, 'dopant sweep'))
    calc_dopant_sweep(conditions, defects, masters)
else
    error('select a calc command')
end

%%%%%%%%%%%%% End of Main %%%%%%%%%%%%%%%

%%%%%%%%%%%%% Subroutines %%%%%%%%%%%%%%%

function calc_single(conditions, defects, masters)
[equilib_dark_sol] = defect_equilibrium_dark(conditions, defects);
    KROGER_Save_Full_Equlibrium_Outputs

    [fullquench_dark_sol] = defect_fullquench_dark(equilib_dark_sol, conditions, defects);
    KROGER_Save_Full_Quench_Outputs

    KROGER_Plot_Outputs
    drawnow
    
    if strcmp(masters.auto_save_flag, 'on')
       masters.figs = [h1, h2, h3, h4, h5, h6, h7, h8];
       save_selected_figures(conditions, defects, masters)
    end
end

function calc_N_sweep(conditions, defects, masters)
for abc = 1:masters.num_N_sweep

    if strcmp(masters.N_sweep_flag, 'Na')
        conditions.Na = masters.N_sweep(abc);
    elseif strcmp(masters.N_sweep_flag, 'Nd')
        conditions.Nd = masters.N_sweep(abc);
    else
        error('when doing Na/Nd sweep you must flag either Na or Nd')
    end

    [equilib_dark_sol] = defect_equilibrium_dark(conditions, defects);
    KROGER_Save_Full_Equlibrium_Outputs

    [fullquench_dark_sol] = defect_fullquench_dark(equilib_dark_sol, conditions, defects);
    KROGER_Save_Full_Quench_Outputs

    KROGER_Plot_Outputs
    drawnow
    
    if strcmp(masters.auto_save_flag, 'on')
       masters.figs = [h1, h2, h3, h4, h5, h6, h7, h8];
       save_selected_figures(conditions, defects, masters)
    end

    pause

end
end

function calc_dopant_sweep(conditions, defects, masters)
for abc = 1:masters.num_dopant_sweep

    conditions.fixed_conc_values(masters.dopant_to_sweep_number) = masters.dopant_sweep(abc);

    [equilib_dark_sol] = defect_equilibrium_dark(conditions, defects);
    KROGER_Save_Full_Equlibrium_Outputs

    [fullquench_dark_sol] = defect_fullquench_dark(equilib_dark_sol, conditions, defects);
    KROGER_Save_Full_Quench_Outputs

    KROGER_Plot_Outputs
    drawnow
    
    if strcmp(masters.auto_save_flag, 'on')
       masters.figs = [h1, h2, h3, h4, h5, h6, h7, h8];
       save_selected_figures(conditions, defects, masters)
    end

    pause

end
end

function save_selected_figures(conditions, defects, masters)

dopant_string = "";


for i = 1:length(defects.elementnames)

    val = conditions.fixed_conc_values(i);

    if val > 0
        element = defects.elementnames(i);
        label = strrep(sprintf('%s%.0e', element, val), '+','');
        dopant_string = dopant_string + "_" + label;
    end

end

if conditions.Na > 0
    N_label = strrep(sprintf('Na%.0e', conditions.Na), '+','');
elseif conditions.Nd > 0
    N_label = strrep(sprintf('Nd%.0e', conditions.Nd), '+','');
else
    N_label = "N0";
end


N_folder = N_label;

save_folder = fullfile(masters.save_root, N_folder);

if ~exist(save_folder, 'dir')
    mkdir(save_folder)
end

for a = masters.auto_save

    filename = sprintf('fig%d_%s_%s%s.png', a, conditions.mu_conditions_flag, N_label, dopant_string);

    filepath = fullfile(save_folder, filename);

    exportgraphics(masters.figs(a), filepath);

end

end

% %% some automatic input validation before proceeding %%%%%%%%%%%%%%%%%%%%%%
% if size(conditions.muT_equilibrium(1,:),2) ~= size(defects.cs_dm,2)
%     error('number of rows in mu must equal number of columns in cs_dm')
% elseif (size(defects.cs_dHo,1) ~= size(defects.cs_charge,1))
%     error('dE and cs_charge must be the same size')
% elseif  ~= size(defects.cs_dm,2) ~= conditions.num_elements
%     error('cs_dm is a matrix with rows for each defect and number of columns that should = number of elements')
% elseif isscalar(conditions.T_fullquench)
%     error('for now, T_quench needs to be a scalar - vectorize it later?')
% else
%     disp('Inputs passed some basic checks on compatibility.  User MUST still check that they did what they really anted to.  Garbage In = Garbagee Out!')
% end


%% save the problem definition in the conditions variable into a mat file for reference
% save(strcat(conditions.save_pname,conditions.conditions_save_fname,'.mat'),'conditions')
save(strcat(conditions.save_pname,conditions.conditions_save_fname,'.mat'))
disp('Saved problem specificaiton conditions to .mat file')
clear conditions_save_fname


%% do the desired calculations - some loops could be built here and resutls captured and saved.  The basic examples are just full equilibrium and full quenching
% Full Equilibrium Calculation
tic
disp('Starting full equilibrium calculation')
[equilib_dark_sol] = defect_equilibrium_dark_with_frozen(conditions, defects);
toc


% Full Quenching Calculation
tic
disp('Starting full quench calculation')
[fullquench_dark_sol] = defect_fullquench_dark(equilib_dark_sol, conditions, defects);
toc
%% this has to come after a full equilibrium calculation

%% compute the stoichiometry and element totals and write it to a file
element_totals_headers = ['T_equilibrium (K)' defects.elementnames];
element_totals_cell = [cellstr(element_totals_headers); num2cell([equilib_dark_sol.T_equilibrium equilib_dark_sol.element_totals])];
writecell(element_totals_cell,strcat(conditions.save_pname,conditions.elements_save_fname),'Delimiter','tab','WriteMode','overwrite');
disp('wrote element total concentrations to text file')


%% group defects into categories - so if there are 5 flavors of VGa you can put them all together for convenience
% in the .mat file, each row is a list of column numbers from the main defect database to sum into a group 

if strcmp(conditions.defect_group_flag,'On')
    temp = load(strcat(conditions.defect_grouping_pname,conditions.defect_grouping_fname));  % this loads 3 variables into a structue "temp'
    conditions.group_together_columns = temp.group_together_columns;
    conditions.group_names = temp.group_names;
    conditions.num_groups = temp.num_groups;
    grouped_defects = zeros(numel(equilib_dark_sol.T_equilibrium), conditions.num_groups);

    for index = 1:conditions.num_groups
        cols_to_sum = conditions.group_together_columns(index, conditions.group_together_columns(index, :) > 0);
        grouped_defects(:, index) = sum(equilib_dark_sol.defects(:, cols_to_sum), 2);
    end

    equilib_grouped_STH = equilib_dark_sol.sth1 + equilib_dark_sol.sth2;
    quenched_grouped_STH = fullquench_dark_sol.sth1 + fullquench_dark_sol.sth2;

    grouping_headers = ['T_equilibrium (K)' 'n equilib (#/cm3)' 'p equilib (#/cm3)' 'STH equilib (#/cm3)' 'n quench (#/cm3)' 'p quench (#/cm3)' 'STH quench (#/cm3)' conditions.group_names'];
    grouping_defects_out = cat(2,equilib_dark_sol.T_equilibrium, equilib_dark_sol.n, equilib_dark_sol.p, equilib_grouped_STH, fullquench_dark_sol.n, fullquench_dark_sol.p,  quenched_grouped_STH, grouped_defects);
    grouping_defects_cell = [cellstr(grouping_headers); num2cell(grouping_defects_out)];

    writecell(grouping_defects_cell,strcat(conditions.save_pname, conditions.grouping_save_fname),'Delimiter','tab','WriteMode','overwrite');  % write results with column headers

    clear temp grouping_save_fname grouping_defects_cell grouping_defects_out grouping_headers grouped_defects index cols_to_sum grouped_STH
end





%%%%%%%%%%% save the equilibrium outputs
all_defects_out = cat(2,equilib_dark_sol.T_equilibrium, equilib_dark_sol.charge_bal_err, equilib_dark_sol.element_bal_err, equilib_dark_sol.tot_bal_err, equilib_dark_sol.mu, equilib_dark_sol.element_totals, conditions.EvT_equilibrium', conditions.EcT_equilibrium', equilib_dark_sol.EFn, equilib_dark_sol.EFp, equilib_dark_sol.Nd, equilib_dark_sol.Na, equilib_dark_sol.n, equilib_dark_sol.p, equilib_dark_sol.sth1, equilib_dark_sol.sth2, equilib_dark_sol.defects);
all_chargestates_out = cat(2,equilib_dark_sol.T_equilibrium, equilib_dark_sol.charge_bal_err, equilib_dark_sol.element_bal_err, equilib_dark_sol.tot_bal_err, equilib_dark_sol.mu, equilib_dark_sol.element_totals, conditions.EvT_equilibrium', conditions.EcT_equilibrium', equilib_dark_sol.EFn, equilib_dark_sol.EFp, equilib_dark_sol.Nd, equilib_dark_sol.Na, equilib_dark_sol.n, equilib_dark_sol.p, equilib_dark_sol.sth1, equilib_dark_sol.sth2, equilib_dark_sol.chargestates);
all_defect_headers = ['T_equilibrium (K)' 'charge_bal_err' 'element_bal_err' 'tot_bal_err' defects.mu_names_with_units defects.elementnames 'Ev(T) (eV)' 'Ec(T) (eV)' 'Efn(eV)' 'Efp (eV)' 'Nd (#/cm3)' 'Na (#/cm3)' 'n (#/cm3)' 'p (#/cm3)' 'sth1 (#/cm3)' 'sth2 (#/cm3)' defects.defect_names'];
all_chargestate_headers = ['T_equilibrium (K)' 'charge_bal_err' 'element_bal_err' 'tot_bal_err' defects.mu_names_with_units defects.elementnames 'Ev(T) (eV)' 'Ec(T) (eV)' 'Efn(eV)' 'Efp (eV)' 'Nd (#/cm3)' 'Na (#/cm3)' 'n (#/cm3)' 'p (#/cm3)' 'sth1 (#/cm3)' 'sth2 (#/cm3)' defects.chargestate_names'];
all_defects_cell = [cellstr(all_defect_headers); num2cell(all_defects_out)];
all_chargestates_cell = [cellstr(all_chargestate_headers); num2cell(all_chargestates_out)];
all_elements_out = cat(2,equilib_dark_sol.T_equilibrium, equilib_dark_sol.element_totals);

% throw out all the ones with insignificant concentrations and save
sig_defects_index = find(max(equilib_dark_sol.defects,[],1)>=conditions.save_min);
sig_chargestates_index = find(max(equilib_dark_sol.chargestates,[],1)>=conditions.save_min);
sig_defects_out = cat(2,equilib_dark_sol.T_equilibrium, equilib_dark_sol.charge_bal_err, equilib_dark_sol.element_bal_err, equilib_dark_sol.tot_bal_err, equilib_dark_sol.mu, conditions.EvT_equilibrium', conditions.EcT_equilibrium', equilib_dark_sol.EFn, equilib_dark_sol.EFp, equilib_dark_sol.Nd, equilib_dark_sol.Na, equilib_dark_sol.n, equilib_dark_sol.p, equilib_dark_sol.sth1, equilib_dark_sol.sth2, equilib_dark_sol.defects(:,sig_defects_index));
sig_chargestates_out = cat(2,equilib_dark_sol.T_equilibrium, equilib_dark_sol.charge_bal_err, equilib_dark_sol.element_bal_err, equilib_dark_sol.tot_bal_err, equilib_dark_sol.mu, conditions.EvT_equilibrium', conditions.EcT_equilibrium', equilib_dark_sol.EFn, equilib_dark_sol.EFp, equilib_dark_sol.Nd, equilib_dark_sol.Na, equilib_dark_sol.n, equilib_dark_sol.p, equilib_dark_sol.sth1, equilib_dark_sol.sth2, equilib_dark_sol.chargestates(:,sig_chargestates_index));
sig_defect_headers = ['T_equilibrium (K)' 'charge_bal_err' 'element_bal_err' 'tot_bal_err' defects.mu_names_with_units 'Ev(T) (eV)' 'Ec(T) (eV)' 'Efn(eV)' 'Efp (eV)' 'Nd (#/cm3)' 'Na (#/cm3)' 'n (#/cm3)' 'p (#/cm3)' 'sth1 (#/cm3)' 'sth2 (#/cm3)' defects.defect_names(sig_defects_index)'];
sig_chargestate_headers = ['T_equilibrium (K)' 'charge_bal_err' 'element_bal_err' 'tot_bal_err' defects.mu_names_with_units 'Ev(T) (eV)' 'Ec(T) (eV)' 'Efn(eV)' 'Efp (eV)' 'Nd (#/cm3)' 'Na (#/cm3)' 'n (#/cm3)' 'p (#/cm3)' 'sth1 (#/cm3)' 'sth2 (#/cm3)' defects.chargestate_names(sig_chargestates_index)'];
sig_defects_cell = [cellstr(sig_defect_headers); num2cell(sig_defects_out)];
sig_chargestates_cell = [cellstr(sig_chargestate_headers); num2cell(sig_chargestates_out)];

% write output tables to a tab delim text file
if strcmp(conditions.save_files_flag,'All')
    writecell(all_defects_cell,strcat(conditions.save_pname, conditions.equilib_save_fname,'_all_defects'),'Delimiter','tab','WriteMode','overwrite');  % write results with column headers
    writecell(all_chargestates_cell,strcat(conditions.save_pname, conditions.equilib_save_fname,'_all_chargestates'),'Delimiter','tab','WriteMode','overwrite');
elseif strcmp(conditions.save_files_flag,'Sig_Only')
    writecell(sig_defects_cell,strcat(conditions.save_pname, conditions.equilib_save_fname,'_sig_defects'),'Delimiter','tab','WriteMode','overwrite');
    writecell(sig_chargestates_cell,strcat(conditions.save_pname, conditions.equilib_save_fname,'_sig_chargestates'),'Delimiter','tab','WriteMode','overwrite');
end
disp('Wrote equilibrium dark solution to text files')

% write the whole solution structure to a mat file
save(strcat(conditions.save_pname,conditions.equilib_save_fname))
disp('Wrote equilibrium dark solution to .mat file')

% clean up memory
clear equilib_save_fname all_defects_out all_chargestates_out all_defects_headers all_chargestates_headers all_defects_cell all_chargstates_cell sig_defects_index sig_chargestates_index sig_defects_out sig_chargestates_out sig_defects_headers sig_chargestates_headers sig_defects_cell sig_chargestates_cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% end full equilibrium calc %%%
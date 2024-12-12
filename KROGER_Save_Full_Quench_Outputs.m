%% Task Trange 2: quenching from the annealing temperatures from example 1  %%%%%%%%%%%%%%%

%%% note this assumes equilib_dark_sol is still in memory - could also
%%% add a line here to load a prior solution from a mat file.


% these variables don't exist for quenching so set them to zero to keep the
% file format the same with the full equilibrium calc
fullquench_dark_sol.element_bal_err = zeros(size(fullquench_dark_sol.charge_bal_err));
fullquench_dark_sol.tot_bal_err = zeros(size(fullquench_dark_sol.charge_bal_err));
fullquench_dark_sol.mu = zeros(size(equilib_dark_sol.mu));
EvT_fullquench_vec = conditions.EvT_fullquench*ones(size(fullquench_dark_sol.charge_bal_err));
EcT_fullquench_vec = conditions.EcT_fullquench*ones(size(fullquench_dark_sol.charge_bal_err));

%%%%%%%%%%% save the outputs for quenching.  Yes some of these are the
%%%%%%%%%%% same as computed for the equilibrium solution but for sake
%%%%%%%%%%% of debugging/proofreading were jsut doing it here again
all_defects_out = cat(2,fullquench_dark_sol.T_equilibrium, fullquench_dark_sol.charge_bal_err, fullquench_dark_sol.element_bal_err, fullquench_dark_sol.tot_bal_err, fullquench_dark_sol.mu, conditions.EvT_equilibrium', conditions.EcT_equilibrium', fullquench_dark_sol.EFn, fullquench_dark_sol.EFp, fullquench_dark_sol.Nd, fullquench_dark_sol.Na, fullquench_dark_sol.n, fullquench_dark_sol.p, fullquench_dark_sol.sth1, fullquench_dark_sol.sth2, fullquench_dark_sol.defects);
all_chargestates_out = cat(2,fullquench_dark_sol.T_equilibrium, fullquench_dark_sol.charge_bal_err, fullquench_dark_sol.element_bal_err, fullquench_dark_sol.tot_bal_err, fullquench_dark_sol.mu, conditions.EvT_equilibrium', conditions.EcT_equilibrium', fullquench_dark_sol.EFn, fullquench_dark_sol.EFp, fullquench_dark_sol.Nd, fullquench_dark_sol.Na, fullquench_dark_sol.n, fullquench_dark_sol.p, fullquench_dark_sol.sth1, fullquench_dark_sol.sth2, fullquench_dark_sol.chargestates);
all_defect_headers = ['T_equilibrium (K)' 'charge_bal_err' 'element_bal_err' 'tot_bal_err' defects.mu_names_with_units 'Ev(T) (eV)' 'Ec(T) (eV)' 'Efn(eV)' 'Efp (eV)' 'Nd (#/cm3)' 'Na (#/cm3)' 'n (#/cm3)' 'p (#/cm3)' 'sth1 (#/cm3)' 'sth2 (#/cm3)' defects.defect_names'];
all_chargestate_headers = ['T_equilibrium (K)' 'charge_bal_err' 'element_bal_err' 'tot_bal_err' defects.mu_names_with_units 'Ev(T) (eV)' 'Ec(T) (eV)' 'Efn(eV)' 'Efp (eV)' 'Nd (#/cm3)' 'Na (#/cm3)' 'n (#/cm3)' 'p (#/cm3)' 'sth1 (#/cm3)' 'sth2 (#/cm3)' defects.chargestate_names'];
all_defects_cell = [cellstr(all_defect_headers); num2cell(all_defects_out)];
all_chargestates_cell = [cellstr(all_chargestate_headers); num2cell(all_chargestates_out)];

% throw out all the ones with insignificant concentrations and save
sig_defects_index = find(max(fullquench_dark_sol.defects,[],1)>=conditions.save_min);
sig_chargestates_index = find(max(fullquench_dark_sol.chargestates,[],1)>=conditions.save_min);
sig_defects_out = cat(2,fullquench_dark_sol.T_equilibrium, fullquench_dark_sol.charge_bal_err, fullquench_dark_sol.element_bal_err, fullquench_dark_sol.tot_bal_err, fullquench_dark_sol.mu, conditions.EvT_equilibrium', conditions.EcT_equilibrium', fullquench_dark_sol.EFn, fullquench_dark_sol.EFp, fullquench_dark_sol.Nd, fullquench_dark_sol.Na, fullquench_dark_sol.n, fullquench_dark_sol.p, fullquench_dark_sol.sth1, fullquench_dark_sol.sth2, fullquench_dark_sol.defects(:,sig_defects_index));
sig_chargestates_out = cat(2,fullquench_dark_sol.T_equilibrium, fullquench_dark_sol.charge_bal_err, fullquench_dark_sol.element_bal_err, fullquench_dark_sol.tot_bal_err, fullquench_dark_sol.mu, conditions.EvT_equilibrium', conditions.EcT_equilibrium', fullquench_dark_sol.EFn, fullquench_dark_sol.EFp, fullquench_dark_sol.Nd, fullquench_dark_sol.Na, fullquench_dark_sol.n, fullquench_dark_sol.p, fullquench_dark_sol.sth1, fullquench_dark_sol.sth2, fullquench_dark_sol.chargestates(:,sig_chargestates_index));
sig_defects_headers = ['T_equilibrium (K)' 'charge_bal_err' 'element_bal_err' 'tot_bal_err' defects.mu_names_with_units 'Ev(T) (eV)' 'Ec(T) (eV)' 'Efn(eV)' 'Efp (eV)' 'Nd (#/cm3)' 'Na (#/cm3)' 'n (#/cm3)' 'p (#/cm3)' 'sth1 (#/cm3)' 'sth2 (#/cm3)' defects.defect_names(sig_defects_index)'];
sig_chargestate_headers = ['T_equilibrium (K)' 'charge_bal_err' 'element_bal_err' 'tot_bal_err' defects.mu_names_with_units 'Ev(T) (eV)' 'Ec(T) (eV)' 'Efn(eV)' 'Efp (eV)' 'Nd (#/cm3)' 'Na (#/cm3)' 'n (#/cm3)' 'p (#/cm3)' 'sth1 (#/cm3)' 'sth2 (#/cm3)' defects.chargestate_names(sig_chargestates_index)'];
sig_defects_cell = [cellstr(sig_defect_headers); num2cell(sig_defects_out)];
sig_chargestates_cell = [cellstr(sig_chargestate_headers); num2cell(sig_chargestates_out)];

% write to a tab delim text file
if strcmp(conditions.save_files_flag,'All')
    writecell(all_defects_cell,strcat(conditions.save_pname,conditions.fullquench_save_fname,'_all_defects'),'Delimiter','tab','WriteMode','overwrite');  % write results with column headers
    writecell(all_chargestates_cell,strcat(conditions.save_pname,conditions.fullquench_save_fname,'_all_chargestates'),'Delimiter','tab','WriteMode','overwrite');
elseif strcmp(conditions.save_files_flag,'Sig_Only')
    writecell(sig_defects_cell,strcat(conditions.save_pname,conditions.fullquench_save_fname,'_sig_defects'),'Delimiter','tab','WriteMode','overwrite');
    writecell(sig_chargestates_cell,strcat(conditions.save_pname,conditions.fullquench_save_fname,'_sig_chargestates'),'Delimiter','tab','WriteMode','overwrite');
end
disp('Wrote full quench dark solution as text files')

% write the whole solution structure to a mat file
save(strcat(conditions.save_pname,conditions.fullquench_save_fname))
disp('Wrote full quench dark solution to .mat file')

% clean up memory
clear fullquench_save_fname all_defects_out all_chargestates_out all_defects_headers all_chargestates_headers all_defects_cell all_chargstates_cell sig_defects_index sig_chargestates_index sig_defects_out sig_chargestates_out sig_defects_headers sig_chargestates_headers sig_defects_cell sig_chargestates_cell  EvT_fullquench_vec EcT_fullquench_vec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% end quenched calc %%%%%%%%%%%%%%%%%
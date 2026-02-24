kBT = 0.050;
cs_concs  = 1e16*ones(1,52);
conditions.num_sites = 5e22*ones(1,4);


% [site_occ_fracs] = site_occupation_fracs(conditions, defects, cs_concs);



[Do_maxD, Emig_maxD, D_maxD, D_tot] = calc_Diff_consts(cs_concs, conditions, defects, kBT);  
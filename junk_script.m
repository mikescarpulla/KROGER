kBT = 0.050;
cs_concs  = 1e10*ones(1,3);
conditions.num_sites = [1e22 2e22];
defects.cs_num_each_site = [1 0; 0 1; 1 1];

[site_occ_fracs] = site_occupation_fracs(conditions, defects, cs_concs);


[Do_maxD, Emig_maxD, D_maxD, D_tot] = calc_Diff_consts(cs_concs, conditions, defects, kBT);  
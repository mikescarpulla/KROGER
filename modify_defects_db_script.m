cs_idx = [1:8 12 13];

defects.cs_names = defects.cs_names(cs_idx);
defects.cs_ID = defects.cs_ID(cs_idx);
defects.cs_q = defects.cs_q(cs_idx); 
defects.cs_Eform = defects.cs_Eform(cs_idx); 
defects.cs_Do = defects.cs_Do(cs_idx); 
defects.cs_Emigration = defects.cs_Emigration(cs_idx); 
defects.cs_degen_factor_config = defects.cs_degen_factor_config(cs_idx);
defects.cs_degen_factor_elec = defects.cs_degen_factor_elec(cs_idx);
defects.cs_site_prefactor = defects.cs_site_prefactor(cs_idx);
defects.cs_tot_prefactor = defects.cs_tot_prefactor(cs_idx);

defects.cs_dm = defects.cs_dm(cs_idx,:); 
defects.cs_num_each_site = defects.cs_num_each_site(cs_idx,:);

def_idx = [1 2 3 5];
defects.defect_names = defects.defect_names(def_idx);
defects.cs_indices_lo = defects.cs_indices_lo(def_idx);
defects.cs_indices_hi = defects.cs_indices_hi(def_idx);
defects.defect_num_cs_per = defects.defect_num_cs_per(def_idx);
           
defects.num_defects = 4; 
defects.num_cs= 10;
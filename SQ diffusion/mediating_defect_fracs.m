function [cs_fracs] = mediating_defect_fracs(defects, cs_concs)

num_D_modes = size(defects.cs_D_mediated_by_cs,2); % in the defect database, how many diffusion modes are allowed for each chargestate?  If one defect has 2 modes, then the whole set of chargestates will have 2 columns for 2 modes, but most will be zero

% check that the sizes of all the variables are right for the case of
% having multiple diffusion modes
if size(defects.Doo,2)~=num_D_modes && size(defects.Emig,2)~=num_D_modes && size(defects.cs_D_mediated_by_cs,2)~=num_D_modes && size(defects.cs_D_mediated_by_site,2)~=num_D_modes
    error('For each diffusion mode, Doo, E_mig, D_mediated_by_cs, and D_mediated_by_site must all have one entry per chargestate as rows and the same number of columns as eachother')
end

cs_fracs = (cs_concs ./ defects.cs_prefactor) * ones(1,num_D_modes);   % divide the conc of each chargestate by its max possible number, then copy this columvector to as many columns as there are diffusion modes.  

% The size of this output variable will be cs x num_diffusion_modes 

end

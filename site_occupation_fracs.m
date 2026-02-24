function [site_occ_fracs] = site_occupation_fracs(conditions, defects, cs_concs)
% This function computes x_jq = N_jq / N_ref for all chargestates in a calculation.  
% conditions.num_sites: this is a vector like [5e22 5e22...]
% cs_num_each_site - this is the matrix of integers telling how many of
% each site is needed to make each chargestate.  
% cs_concs = vector of concentrations 

if ~isvector(cs_concs)
    error('cs_concs needs to be a vector')
end

cs_concs = cs_concs(:)';  % turn this into a row vector which ever one it comes in as
% num_sites_in_crystal = numel(conditions.num_sites);
tot_sites_occupied = cs_concs * defects.cs_num_each_site;   % num each site is a num_cs x num_sites matrix, while cs_concs is a row vector.  


% if the lin alg method doesnt work, here is a for loop 
% tot_sites_occupied = zeros(1,num_sites_in_crystal);  
% for i = 1:num_sites_in_crystal
%     tot_sites_occupied(i) = dot(cs_concs, defects.cs_num_each_site(:,i));
% end



site_occ_fracs = tot_sites_occupied ./ conditions.num_sites;

end

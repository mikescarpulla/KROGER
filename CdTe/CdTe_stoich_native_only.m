function [ratio_stoich, element_totals] = CdTe_stoich_native_only(solution, conditions, defects)


num_T = size(solution.cs, 1);
num_Cd = zeros(num_T,1);
num_Te = zeros(num_T,1);
Cd_Te_stoich = zeros(num_T,2);


for j = 1:num_T  % loops over all the T's
    num_Cd(j) = conditions.num_sites(1) + dot(defects.cs_dm(:,1),solution.cs(j,:));  % dm holds the number of atoms added (+) or removed (-) relative to full latice.  So start with full sublattice of Cd and add the change to get the final #
    num_Te(j) = conditions.num_sites(2) + dot(defects.cs_dm(:,2),solution.cs(j,:));
end

Cd_Te_stoich = [num_Cd num_Te]/1.47e22;

X_Te = num_Te ./ (num_Cd + num_Te);

ratio_stoich = cat(2, X_Te, Cd_Te_stoich);

element_totals = [num_Cd num_Te];



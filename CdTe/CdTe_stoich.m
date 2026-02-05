function [ratio_stoich, element_totals] = CdTe_stoich(solution, conditions, defects)


num_entries = size(solution.chargestates, 1);
num_Cd = zeros(num_entries,1);
num_Te = zeros(num_entries,1);
num_N = zeros(num_entries,1);
num_P = zeros(num_entries,1);
num_As = zeros(num_entries,1);
num_Sb = zeros(num_entries,1);
num_Cu = zeros(num_entries,1);
num_Cl = zeros(num_entries,1);
num_O = zeros(num_entries,1);
Cd_Te_stoich = zeros(num_entries,2);


for j = 1:num_entries  % loops over all the T's
    num_Cd(j) = conditions.num_sites(1) + dot(defects.cs_dm(:,1),solution.chargestates(j,:));
    num_Te(j) = conditions.num_sites(2) + dot(defects.cs_dm(:,2),solution.chargestates(j,:));
    num_N(j) = dot(defects.cs_dm(:,3),solution.chargestates(j,:));
    num_P(j) = dot(defects.cs_dm(:,4),solution.chargestates(j,:));
    num_As(j) = dot(defects.cs_dm(:,5),solution.chargestates(j,:));
    num_Sb(j) = dot(defects.cs_dm(:,6),solution.chargestates(j,:));
    num_Cu(j) = dot(defects.cs_dm(:,7),solution.chargestates(j,:)); 
    num_Cl(j) = dot(defects.cs_dm(:,8),solution.chargestates(j,:));
    num_O(j) = dot(defects.cs_dm(:,9),solution.chargestates(j,:)); 


    Cd_Te_stoich(j,:) = [num_Cd(j) num_Te(j)]/1.47e22;
end

ratio = Cd_Te_stoich(:,2)./(Cd_Te_stoich(:,1) + Cd_Te_stoich(:,2));

ratio_stoich = cat(2,ratio,Cd_Te_stoich);

element_totals = [num_Cd num_Te num_N num_P num_As num_Sb num_Cu num_Cl num_O];



function [ratio_stoich, element_totals] = GaN_stoich(solution, conditions, defects)


num_entries = size(solution.chargestates, 1);
num_Ga = zeros(num_entries,1);
num_N = zeros(num_entries,1);
num_C = zeros(num_entries,1);
num_Si = zeros(num_entries,1);
num_Mg = zeros(num_entries,1);
num_Ge = zeros(num_entries,1);
num_O = zeros(num_entries,1);
num_Be = zeros(num_entries,1);
num_Zn = zeros(num_entries,1);
num_Ca = zeros(num_entries,1);
num_Cd = zeros(num_entries,1);
num_H = zeros(num_entries,1);
Ga_N_stoich = zeros(num_entries,2);


for j = 1:num_entries  % loops over all the T's
    num_Ga(j) = conditions.num_sites(1) + dot(defects.cs_dm(:,1),solution.chargestates(j,:));
    num_N(j) = conditions.num_sites(2) + dot(defects.cs_dm(:,2),solution.chargestates(j,:));
    num_C(j) = dot(defects.cs_dm(:,3),solution.chargestates(j,:));
    num_Si(j) = dot(defects.cs_dm(:,4),solution.chargestates(j,:));
    num_Mg(j) = dot(defects.cs_dm(:,5),solution.chargestates(j,:));
    num_Ge(j) = dot(defects.cs_dm(:,6),solution.chargestates(j,:));
    num_O(j) = dot(defects.cs_dm(:,7),solution.chargestates(j,:));
    num_Be(j) = dot(defects.cs_dm(:,8),solution.chargestates(j,:));
    num_Zn(j) = dot(defects.cs_dm(:,9),solution.chargestates(j,:));
    num_Ca(j) = dot(defects.cs_dm(:,10),solution.chargestates(j,:));
    num_Cd(j) = dot(defects.cs_dm(:,11),solution.chargestates(j,:));
    num_H(j) = dot(defects.cs_dm(:,12),solution.chargestates(j,:));


    Ga_N_stoich(j,:) = [num_Ga(j) num_N(j)]/4.37e22;
end

ratio = Ga_N_stoich(:,2)./(Ga_N_stoich(:,1) + Ga_N_stoich(:,2));

ratio_stoich = cat(2,ratio,Ga_N_stoich);

element_totals = [num_Ga num_N num_C num_Si num_Mg num_Ge num_O num_Be num_Zn num_Ca num_Cd num_H];



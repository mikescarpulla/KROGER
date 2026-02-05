function [ratio_stoich, element_totals] = Ga2O3_stoich(solution, conditions, defects)
% function [Ga_O_stoich, num_Ga, num_O, num_Si, num_H, num_Fe, num_Sn, num_Cr, num_Ti, num_Ir, num_Mg, num_Ca, num_Zn, num_Co, num_Zr, num_Hf, num_Ta, num_Ge] =  Ga2O3_stoich(solution, conditions, defects)


num_entries = size(solution.chargestates, 1);
num_Ga = zeros(num_entries,1);
num_O = zeros(num_entries,1);
num_Si = zeros(num_entries,1);
num_H = zeros(num_entries,1);
num_Fe = zeros(num_entries,1);
num_Sn = zeros(num_entries,1);
num_Cr = zeros(num_entries,1);
num_Ti = zeros(num_entries,1);
num_Ir = zeros(num_entries,1);
num_Mg = zeros(num_entries,1);
num_Ca = zeros(num_entries,1);
num_Zn = zeros(num_entries,1);
num_Co = zeros(num_entries,1);
num_Zr = zeros(num_entries,1);
num_Hf = zeros(num_entries,1);
num_Ta = zeros(num_entries,1);
num_Ge = zeros(num_entries,1);
num_Pt = zeros(num_entries,1);
num_Rh = zeros(num_entries,1);
Ga_O_stoich = zeros(num_entries,2);



for j = 1:num_entries  % loops over all the T's
    num_Ga(j) = conditions.num_sites(1) + conditions.num_sites(2) + dot(defects.cs_dm(:,1),solution.chargestates(j,:));
    num_O(j) = conditions.num_sites(3) + conditions.num_sites(4) + conditions.num_sites(5) + dot(defects.cs_dm(:,2),solution.chargestates(j,:));
    num_Si(j) = dot(defects.cs_dm(:,3),solution.chargestates(j,:));
    num_H(j) = dot(defects.cs_dm(:,4),solution.chargestates(j,:));
    num_Fe(j) = dot(defects.cs_dm(:,5),solution.chargestates(j,:));
    num_Sn(j) = dot(defects.cs_dm(:,6),solution.chargestates(j,:));
    num_Cr(j) = dot(defects.cs_dm(:,7),solution.chargestates(j,:));     
    num_Ti(j) = dot(defects.cs_dm(:,8),solution.chargestates(j,:));
    num_Ir(j) = dot(defects.cs_dm(:,9),solution.chargestates(j,:));
    num_Mg(j) = dot(defects.cs_dm(:,10),solution.chargestates(j,:));
    num_Ca(j) = dot(defects.cs_dm(:,11),solution.chargestates(j,:));
    num_Zn(j) = dot(defects.cs_dm(:,12),solution.chargestates(j,:));
    num_Co(j) = dot(defects.cs_dm(:,13),solution.chargestates(j,:));
    num_Zr(j) = dot(defects.cs_dm(:,14),solution.chargestates(j,:));
    num_Hf(j) = dot(defects.cs_dm(:,15),solution.chargestates(j,:));
    num_Ta(j) = dot(defects.cs_dm(:,16),solution.chargestates(j,:));
    num_Ge(j) = dot(defects.cs_dm(:,17),solution.chargestates(j,:));
    num_Pt(j) = dot(defects.cs_dm(:,18),solution.chargestates(j,:));
    num_Rh(j) = dot(defects.cs_dm(:,19),solution.chargestates(j,:));
    Ga_O_stoich(j,:) = [num_Ga(j) num_O(j)]/1.91e22;
end

ratio = Ga_O_stoich(:,2)./(Ga_O_stoich(:,1) + Ga_O_stoich(:,2));

ratio_stoich = cat(2, ratio, Ga_O_stoich);

element_totals = [num_Ga num_O num_Si num_H num_Fe num_Sn num_Cr num_Ti num_Ir num_Mg num_Ca num_Zn num_Co num_Zr num_Hf num_Ta num_Ge num_Pt num_Rh];



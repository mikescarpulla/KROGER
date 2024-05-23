function [Ga2O3_stable_mask, Ga2O3_decomposes_mask, instability_flags] = Ga2O3_stability(T,Ptot,pGa_vap,pO2,pGa2O,pGaO,P_units) 

instability_flags = zeros(1,4);   % the flag will be 1 for each reaction if it makes Ga2O3 decompose

%%%%%% Check if the specified conditions make Ga2O3
%%%%%% stable or unstable to each reaction.  

% Rxn1 = decomposition of Ga2O3 into liq/solid Ga plus O2 gas.  Ga203(s) = 2Ga(ls) + 3/2 O2(g)
rxn1_dG0 = 2*Ga_ls_G0(T,1,'atm') + 3/2*O2_G0(T,1,'atm') - Ga2O3_G0(T,1,'atm');
rxn1_K0 = exp(rxn1_dG0./(kB_eV*T));
rxn1_K = pO2.^(3/2);
rxn1_unstable_mask = rxn1_K < rxn1_K0;
rxn1_stable_mask = rxn1_K >= rxn1_K0;
if sum(rxn1_unstable_mask)>0
    instability_flags(1) = 1;
end

% Rxn2 = decomposition of Ga2O3 into vapor Ga plus O2 gas.  Ga203(s) = 2Ga(vap) + 3/2 O2(g)
rxn2_dG0 = 2*Ga_vapor_G0(T,1,'atm') + 3/2*O2_G0(T,1,'atm') - Ga2O3_G0(T,1,'atm');
rxn2_K0 = exp(rxn2_dG0./(kB_eV*T));
rxn2_K = pO2.^(3/2) * pGa_vap.^2;
rxn2_unstable_mask = rxn2_K < rxn2_K0;
rxn2_stable_mask = rxn2_K >= rxn1_K0;
if sum(rxn2_unstable_mask)>0
    instability_flags(2) = 1;
end

% Rxn3 = decomposition of Ga2O3 into Ga2O plus O2 gas.  Ga203(s) = 2Ga(vap) + 3/2 O2(g)
rxn3_dG0 = Ga2O_G0(T,1,'atm') + O2_G0(T,1,'atm') - Ga2O3_G0(T,1,'atm');
rxn3_K0 = exp(rxn3_dG0./(kB_eV*T));
rxn3_K = pO2 * pGa2O;
rxn3_unstable_mask = rxn3_K < rxn3_K0;
rxn3_stable_mask = rxn3_K >= rxn3_K0;
if sum(rxn3_unstable_mask)>0
    instability_flags(3) = 1;
end

% Rxn4 = decomposition of Ga2O3 into GaO plus O2 gas.  2Ga203(s) = 4Ga0(vap) + O2(g)
rxn4_dG0 = 4*GaO_G0(T,1,'atm') + O2_G0(T,1,'atm') - 2*Ga2O3_G0(T,1,'atm');
rxn4_K0 = exp(rxn4_dG0./(kB_eV*T));
rxn4_K = pO2 * pGaO.^4;
rxn4_unstable_mask = rxn4_K < rxn4_K0;
rxn4_stable_mask = rxn4_K >= rxn4_K0;
if sum(rxn4_unstable_mask)>0
    instability_flags(4) = 1;
end

% Rxn5 = Ga(l/s) + Ga2O3 decomposes into Ga2O plus O2 gas.  2Ga(l/s) + Ga203(s) = 2Ga20(vap) + 1/2 O2(g)
rxn5_dG0 = 2*Ga2O_G0(T,1,'atm') + 1/2*O2_G0(T,1,'atm') - Ga2O3_G0(T,1,'atm') - 2*Ga_ls_G0(T,1,'atm');
rxn5_K0 = exp(rxn5_dG0./(kB_eV*T));
rxn5_K = pO2.^0.5 * pGa2O.^2;
rxn5_unstable_mask = rxn5_K < rxn5_K0;
rxn5_stable_mask = rxn5_K >= rxn5_K0;
if sum(rxn5_unstable_mask)>0
    instability_flags(4) = 1;
end




Ga2O3_decomposes_mask = zeros(size(T)) < (rxn1_stable_mask + rxn2_stable_mask + rxn3_stable_mask + rxn4_stable_mask + rxn5_stable_mask);   % create the overall mask for whether Ga2O3 will decompose by any means at any of the T's given.  If there is a 1, it means Ga2O3 will decompse at that temperature.  
Ga2O3_stable_mask = ones(size(Ga2O3_decomposes_mask)) - Ga2O3_decomposes_mask;   % the stability mask is 1- the decomposes one.  









    %
    %         instability_mask = (mu_O<mu_O_Garich) .* (mu_O<mu_O_Ga2Orich) ;   %ones and zeros for unstable conditions
    % %         instability_mask = find((mu_O<mu_O_Garich) .* (mu_O>mu_O_Orich) .* (mu_O<mu_O_Ga2Orich)) ;  % use this line if yo uwant to also detect when mu_O is too high - that could only happen if pO2 is larger than Ptot (impossible) or if some enegetic species like a plasma is used
    %         stability_mask = ones(size(instability_mask))-instability_mask;
    %         instability_index = find(instability_mask);
    %         stability_index = find(stability_mask);
    %
    %         if sum(instability_mask)>0
    %             for i = 1:10
    %             disp('WARNING!!!!  calculation done with constant pO2.  Ga2O3 is be unstable for some temperatures becasue pO2 goes outside of Ga-rich and O-rich bounds!!!')
    %             disp(strcat('with the given pO2, Ga2O3 is unstable for_',num2str(conditions.T_equilibrium(instability_index)),'_K'))
    %             disp('T_equilibrium will be truncated to only show the valid T range for the given pO2  !!!')
    %             end
    %         end
    %
    %         conditions.T_equilibrium = conditions.T_equilibrium(stability_index);   %% truncate the T list to only show the ones for which Ga2O3 is stable
    %         % truncate all the other things to the right length
    %         conditions.EgT_equilibrium = conditions.EgT_equilibrium(stability_index);
    %         conditions.EcT_equilibrium = conditions.EcT_equilibrium(stability_index);
    %         conditions.EvT_equilibrium = conditions.EgT_equilibrium(stability_index);
    %         conditions.NcT_equilibrium = conditions.NcT_equilibrium(stability_index);
    %         conditions.NvT_equilibrium = conditions.NvT_equilibrium(stability_index);
    %         conditions.mu_Ga2O3 = conditions.mu_Ga2O3(stability_index);
    %         conditions.muT_equilibrium = conditions.muT_equilibrium(stability_index,:);
    %
    %

    %
    %     elseif strcmp(conditions.mu_conditions_flag,'O-rich') || strcmp(conditions.mu_conditions_flag,'Ga-rich') || strcmp(conditions.mu_conditions_flag,'Ga2O')   %% O and Ga
    %         mu_O = 0;
    %         mu_Ga = 0;  %set these two as zero to make the variable the right size, then reset the values from ZA functions
    %         conditions.muT_equilibrium = ones(size(conditions.T_equilibrium,2),1)*[mu_Ga mu_O mu_Si mu_H mu_Fe mu_Sn];  % initializes all the chem potentials
    %         [conditions.muT_equilibrium(:,1), conditions.muT_equilibrium(:,2)] = Ga2O3_chemical_potentials(conditions.T_equilibrium,conditions.P,conditions.P_units,conditions.mu_conditions_flag);   % muGa, muO are set, all others in columns 3-5 for now are set to zero
    %         % elseif strcmp(conditions.T_dep_mu_flag,'manual')
    %         %     conditions.mu = -1*ones(size(conditions.T_equilibrium,2),conditions.num_elements);  % here one could enter a desired table of mu vs T

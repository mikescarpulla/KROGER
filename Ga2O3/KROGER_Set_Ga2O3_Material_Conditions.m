%% Set Semiconductor Bandstructure and Lattice sites %%%%%%%%%%%%%%%%%%%%%%
% this has to come after the temperature vector is set
%%% Ga2O3 %%

conditions.vibent_T0 = 870; % This is the characteristic To for the average phonon mode in Ga2O3, as determined from the wo that makes the Debeye temperature come out to about 740 K based on the Cv(T) data found online.  To where x = hbar*w0/kBT = T0/T, wherein w0 is the mean phonon energy from the Debeye temperature determined from Cv(T).  

% when T-independent values are used
conditions.TRef = 300;
conditions.EgRef = 4.8;   % this is for constant Eg(T)
conditions.EvRef = 0;
conditions.EcRef = conditions.EgRef;
conditions.NcRef = 1e19;
conditions.NvRef = 5e20;


% T-dependent Eg values
% conditions.Eg0 = 5.1;
conditions.Eg0 = 5.0;
% conditions.Eg0 = 4.9;   % bandgap at 0K in eV.  Using Adrian's number from STEM EELS at PSU - 4.8 eV at 300 K.  So about 4.9 eV at 0K using
conditions.varshini_a = 0.00105248259206626;  % can implement any model you want for Eg(T) as long as it gives a value for all T_equilibrium
conditions.varshini_b = 676.975385507958;

conditions.EcT_fraction = 0;
% conditions.EcT_fraction = 0.25;
% conditions.EcT_fraction = 0.50;     % what fraction of Eg(T)  happens in the CB?
% conditions.EcT_fraction = 0.75;
% conditions.EcT_fraction = 1;
conditions.EvT_fraction = 1-conditions.EcT_fraction;

% these take the equilibrium T values from conditions and generate the
% needed T-dependent values.  The fraction of Eg change caused by VB vs CB
% can be modified above.
conditions.NcT_equilibrium = conditions.NcRef*(conditions.T_equilibrium/conditions.TRef).^1.5;
conditions.NvT_equilibrium = conditions.NvRef*(conditions.T_equilibrium/conditions.TRef).^1.5;
delta_EgT_equilibrium = (conditions.varshini_a .* conditions.T_equilibrium.^2)./(conditions.T_equilibrium + conditions.varshini_b);  %delta is a positive number
conditions.EgT_equilibrium = conditions.Eg0 - delta_EgT_equilibrium;
conditions.EcT_equilibrium = conditions.Eg0*ones(size(conditions.EgT_equilibrium)) - conditions.EcT_fraction*delta_EgT_equilibrium;
conditions.EvT_equilibrium = zeros(size(conditions.EgT_equilibrium)) + conditions.EvT_fraction*delta_EgT_equilibrium;

conditions.NcT_fullquench = conditions.NcRef*(conditions.T_fullquench/conditions.TRef).^1.5;
conditions.NvT_fullquench = conditions.NvRef*(conditions.T_fullquench/conditions.TRef).^1.5;
delta_EgT_fullquench = (conditions.varshini_a .* conditions.T_fullquench.^2)./(conditions.T_fullquench + conditions.varshini_b);  %delta is a positive number
conditions.EgT_fullquench = conditions.Eg0 - delta_EgT_fullquench;
conditions.EcT_fullquench = conditions.Eg0 - conditions.EcT_fraction*delta_EgT_fullquench;
conditions.EvT_fullquench = conditions.EvT_fraction*delta_EgT_fullquench;


clear delta_EgT_equilibrium delta_EgT_fullquench   %clean up memory after using these

% site densities in the Ga2O3 lattice (num FU/unit cell / vol for unit cell)
N_Ga1 = 1.91e22;    %site 1
N_Ga2 = 1.91e22;  % site 2
N_O1 = 1.91e22;  %site 3
N_O2 = 1.91e22;   %site 4
N_O3 = 1.91e22;  %site5
N_ia = 1.91e22;
N_ib = 1.91e22;
N_ic = 1.91e22;
N_unspecified_i = 1.91e22;  % interstitials - watch the degeneracy of these
conditions.num_sites = [N_Ga1; N_Ga2; N_O1; N_O2; N_O3; N_ia; N_ib; N_ic; N_unspecified_i];   % note this needs to be a column vector, dont change it to a row vector

% % Calculate the numerical prefactor for each defect from defects.degen_factor_config, defects.degen_factor_elec, defects.cs_num_each_site, and conditions.num_sites

%% defects.cs_num_each_site tells you the number of primitive unit cells needed to form each defect (assuming the sites of the crystal, including distinct interstitials like the ia, ib, ... in b-ga2O3, are counted such that there is one per primitive unit cell
defects.cs_site_prefactor = sort((ones(defects.num_chargestates,1)*conditions.num_sites')./defects.cs_num_each_site,2);      % calculates the max number of each chargestate that could be formed, given the availability of each site type in the crystal (or unit cells, or supercells, whatever the basis for counting).  Then sorts by the rows ascending, so the first column becomes the max number we could form.  Its done with a sort becasue you end up with Inf when you divide by sites you need zero of
defects.cs_site_prefactor = defects.cs_site_prefactor(:,1);  % toss out the other columns keeping only the first one, which should be the limiting one
defects.cs_prefactor = defects.cs_degen_factor_config .* defects.cs_degen_factor_elec .* defects.cs_site_prefactor;

clear N_Ga1 N_Ga2 N_O1 N_O2 N_O3 N_ia N_ib N_ic N_unspecified_i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
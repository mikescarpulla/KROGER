%% compute concentrations of chargestates at the lowest T assuming conc of defects are all frozen at higher T's
% improvements 10/11/2025 by MS:
% simplified repeated blocks of code in chargestate number calculations
% updated definition of the quantum vibent per mode to (n+1)ln(n+1)-n ln(n)
% added dG for defects and chargestates as output variables for both equilibrium and fullquench
% standardized variable names to dummy_X (and got rid of misspelled dumy_X)

function[modified_SQ_dark_sol] = modify_SQ_by_defect_removal(input_SQ_dark_sol, dummy_conditions, dummy_defects,dummy_modification)

% dummy_modification.which_defects = vector list of which ones
% dummy_modification.new_value = vector list of new values


%%%%%  initialize the solution structure.  copy from the equilibrium solution
modified_SQ_dark_sol = input_SQ_dark_sol;

% rmfield clears unnecessary variables from the struct
modified_SQ_dark_sol = rmfield(modified_SQ_dark_sol, {'charge_bal_err','element_bal_err','tot_bal_err'});
modified_SQ_dark_sol = rmfield(modified_SQ_dark_sol, {'n','p','EFn','EFp','chargestates'});

% clear modified_SQ_dark_sol.charge_bal_err modified_SQ_dark_sol.element_bal_err modified_SQ_dark_sol.tot_bal_err
% clear modified_SQ_dark_sol.n modified_SQ_dark_sol.p modified_SQ_dark_sol.EFn modified_SQ_dark_sol.EFp
% clear modified_SQ_dark_sol.chargestates

Nd = input_SQ_dark_sol.Nd(:,:,end);
Na = input_SQ_dark_sol.Na(:,:,end);
clear input_SQ_dark_sol.Na input_SQ_dark_sol.Nd


num_chardist = numel(input_SQ_dark_sol.SQ_chardists);
num_Trates = numel(input_SQ_dark_sol.SQ_Trates);
num_defects = size(input_SQ_dark_sol.defects,4);
num_chargestates = size(input_SQ_dark_sol.chargestates,4);
num_mods = numel(dummy_modification.which_defects);

defects_holder = zeros(num_chardist,num_Trates,num_defects);
chargestates_holder = zeros(num_chardist,num_Trates,num_chargestates);
modified_SQ_dark_sol.n = zeros(num_chardist,num_Trates);
modified_SQ_dark_sol.p = zeros(num_chardist,num_Trates);
modified_SQ_dark_sol.sth1 = zeros(num_chardist,num_Trates);
modified_SQ_dark_sol.sth2 = zeros(num_chardist,num_Trates);
modified_SQ_dark_sol.EF = zeros(num_chardist,num_Trates);
modified_SQ_dark_sol.charge_bal_err = zeros(num_chardist,num_Trates);

% the incoming SQ_sol.defects has indices i,j,l,k where k indexes the defect.  Dont ask why it
% just is.

num_tempcalcs = size(input_SQ_dark_sol.defects,3);   %size of defects along 3rd dim, which stores the

old_defects = input_SQ_dark_sol.defects;
input_SQ_dark_sol.defects = zeros(1,num_defects);
modified_SQ_dark_sol.defects = defects_holder;
modified_SQ_dark_sol.chargestates = chargestates_holder;


for m = 1:num_mods
    old_defects(:,:,:,dummy_modification.which_defects(m)) = dummy_modification.new_value(m);
end


for i_1 = 1:num_chardist      %loop over chardist
    for j = 1:num_Trates    %loop over Trate

        for k = 1:num_defects
            input_SQ_dark_sol.defects(k) =  old_defects(i_1,j,input_SQ_dark_sol.l_value_for_defects_cs(i_1,j),k);
            input_SQ_dark_sol.l_value_for_defects_cs(i_1,j)
        end
        EF_guess = FQ_EF_guess(dummy_conditions);   % get a guess for EF close to minimum using trial and error
        EF_out = fzero(@FQ_charge_bal, EF_guess);
        EF_out = fzero(@FQ_charge_bal, EF_out);  %restart just in case
        modified_SQ_dark_sol.EF(i_1,j) = EF_out;
        [modified_SQ_dark_sol.n(i_1,j), modified_SQ_dark_sol.p(i_1,j), modified_SQ_dark_sol.sth1(i_1,j), modified_SQ_dark_sol.sth2(i_1,j)] = FQ_carrier_concentrations(EF_out);
        modified_SQ_dark_sol.chargestates(i_1,j,:) = FQ_chargestate_concentrations(EF_out);  % compute the defect and carrier concentrations from the EF
        modified_SQ_dark_sol.charge_bal_err(i_1,j) = FQ_charge_bal(EF_out);
    end
end

%%%%%%%%%%%%%% end of main calc  %%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%% subroutines %%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [guess] = FQ_EF_guess(FQEF_conditions)
        % function that uses grid search to get EF close to the charge balance solution.  It will output two values for EF that bracket the solution unless something is strange about the charge_bal vs EF (i.e. its not monotonic   %%%%%%%%%%%%

        EF_int = FQEF_conditions.kBT_fullquench / FQEF_conditions.fullquench_EF_search_step_divisor;   % this guarantees you can't miss the solution whcih should be thus within kB/2 of the guess
        EF_grid = (FQEF_conditions.EvT_fullquench - 5*FQEF_conditions.kBT_fullquench) : EF_int : (ceil(FQEF_conditions.EgT_fullquench/EF_int)*EF_int + 5*FQEF_conditions.kBT_fullquench);  % this makes a grid to check over the range -5kBT to Eg+5kBT.  This is ok since we have Fermi-Dirac stats
        nn = size(EF_grid,2);
        errs = zeros(1,nn);

        for i=1:nn
            guess = EF_grid(i);
            errs(i)=FQ_charge_bal(guess);
        end

        edge_index = find(diff(sign(errs))~=0);   %  this finds the rising/falling edge where error changes sign

        if sum(size(edge_index)==[1 1])==2
            min_index = [edge_index edge_index+1];   % so we are finding the two guesses that bracket the solution one + and one -
            %     elseif sum(size(edge_index)~=[1 1])==2
            %         disp('waring: charge balance error may not be monotonic')
        elseif sum(size(edge_index)==[1 2])==2
            min_index = edge_index;
        else
            error('something strange about charge balance error vs EF - solutions may not be valid')
        end

        guess = EF_grid(min_index);
    end    %%%% end EF_guess



    function [charge_bal] = FQ_charge_bal(EF_dummy)
        %%%% function that computes the net charge
        [n, p, sth1, sth2] = FQ_carrier_concentrations(EF_dummy);
        [N_chargestates_FQ] = FQ_chargestate_concentrations(EF_dummy);
        charge_bal = sum(dummy_defects.cs_charge .* N_chargestates_FQ') + p + sth1 + sth2 - n + dummy_conditions.Nd - dummy_conditions.Na ;
    end  % charge_bal



    function [N_chargestates_FQ] = FQ_chargestate_concentrations(EF_dummy)
        %% Given EF and total number of each defect type, figure out the number of each charge state.

        N_chargestates_FQ = zeros(1,dummy_defects.num_chargestates);% set up arrays
        Z = zeros(1,dummy_defects.num_defects);
        [dG_cs_rel] = dG_chargestates_rel(EF_dummy);
        Boltz_facs = exp(-dG_cs_rel/dummy_conditions.kBT_fullquench);
        for i = 1:dummy_defects.num_defects  % loop over defects (defect.cs_ID) - not over charge states
            indices = dummy_defects.cs_ID == i;  % find the indices of the charge states of the ith defect
            Z(i) = sum(Boltz_facs(indices)); % matrix with Z value for each defect (computed from Boltz factors for each charge state in that defect
            N_chargestates_FQ(indices) = Boltz_facs(indices)/Z(i) * input_SQ_dark_sol.defects(i);     % equilib_dark_sol.defects(i) is a scalar, Z(i) is a scalar
        end
    end  %%%% end chargestate concentrations



    function [dSvib_Q_per_mode_norm] = dSvib_quantum_per_mode(T, T0)
        % T0/T is x = hbar*omega0/kBT.  So compute the T0 for the mean mode
        % from the Debeye temperature.
        % make sure to use T and not kBT!!
        n_ph = 1/(exp(T0/T)-1);
        dSvib_Q_per_mode_norm = (1+n_ph)*log(1+n_ph) - n_ph*log(n_ph)+log(2);
        % dSvib_Q_per_mode_norm = T0/(2*TempK)*coth(T0/(2*TempK)) + log(csch(T0/(2*TempK)));
    end


    function [dG_cs_rel] = dG_chargestates_rel(EF_dummy)
        dG_cs_rel = dummy_defects.cs_dHo + dummy_defects.cs_charge * EF_dummy;  % compute the part without mu and without vibent
        if strcmp(dummy_conditions.vib_ent_flag,'3kB')
            dG_cs_rel = dG_cs_rel - 3*dummy_conditions.kBT_fullquench * sum(dummy_defects.cs_dm,2);  % last term is -TdS.  A vacacny has sum(cs_dm)=-1, and interstitial has +1.  dSvib = 3*kB*sum(cs_dm)*f(T) where f(T) is the classical or quantum function (positive numbers).  dG=dH-TdS = dH - 3kBT*f(T)*sum(dm).
        elseif strcmp(dummy_conditions.vib_ent_flag,'Quantum')
            dG_cs_rel = dG_cs_rel - 3*dummy_conditions.kBT_fullquench * dSvib_quantum_per_mode(dummy_conditions.T_fullquench, dummy_conditions.vibent_T0) * sum(dummy_defects.cs_dm,2);
        elseif strcmp(dummy_conditions.vib_ent_flag,'Off')
            % do nothing here - ignoring dSvib
        else
            error('Vibrational entropy flag must be 3kB, Quantum, or Off')
        end
    end



    function [n, p, sth1, sth2] = FQ_carrier_concentrations(EF_dummy)
        %%%% function to compute n and p.  Input must be EF only (scalar
        %%%% not vector with mu's attached)  %%%%%%%%%%%%
        etaCB = (EF_dummy - dummy_conditions.EcT_fullquench)/dummy_conditions.kBT_fullquench;
        etaVB = (EF_dummy - dummy_conditions.EvT_fullquench)/dummy_conditions.kBT_fullquench;   % this looks wrong (not symmetric compared to CB case) but it is right. Direction of integration and sign on EF-Ev are both swapped.
        etaSTH1 = (EF_dummy - (dummy_conditions.EvT_fullquench + dummy_conditions.E_relax_sth1) )/dummy_conditions.kBT_fullquench;
        etaSTH2 = (EF_dummy - (dummy_conditions.EvT_fullquench + dummy_conditions.E_relax_sth2) )/dummy_conditions.kBT_fullquench;
        if strcmp(dummy_conditions.Boltz_or_FD_flag,'FD')            % use Fermi-Dirac integrals so degenerate conditions handled correctly
            n = n_Fermi_Dirac(etaCB, dummy_conditions.NcT_fullquench);
            p = n_Fermi_Dirac(-etaVB, dummy_conditions.NvT_fullquench);
            sth1 = dummy_conditions.sth_flag * n_Fermi_Dirac(-etaSTH1, dummy_conditions.num_sites(3));
            sth2 = dummy_conditions.sth_flag * n_Fermi_Dirac(-etaSTH2, dummy_conditions.num_sites(4));
        elseif strcmp(dummy_conditions.Boltz_or_FD_flag,'Boltz')              % %     use just Boltzmann approx
            n = dummy_conditions.NcT_fullquench * exp(etaCB);   % these are right (sign swap).  Boltzmann factors should end up <1 when EF is in gap
            p = dummy_conditions.NvT_fullquench * exp(-etaVB);
            sth1 = dummy_conditions.sth_flag * dummy_conditions.num_sites(3) * exp(-etaSTH1);
            sth2 = dummy_conditions.sth_flag * dummy_conditions.num_sites(4) * exp(-etaSTH2);
        else
            error('Boltz_or_FD_flag must be Boltz or FD')
        end
    end



end   %%%%%end main function
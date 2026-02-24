% this function calculates the defect equilibrium Energies entered in eV,
% concentrations entered in and come out in #/cm3 Defect charge state
% formation enthalpy is assumed to be in the form dH = Eform + q*EF don't
% mess with the location of the "end" for the main function - it is placed
% so the subrouties are nested functions and thus can access the workspace
% of the main function - letting them see the conditions, material, and
% defects variables without having those as explicit input arguments.
% 1/14/21  - implemented Fermi Dirac approximations rather than direct
% integration. 12/23/22 - added handling for fixed elements and fixed
% defects, which can then allow generalized quenching calcs 2/7/2023 =
% added in site blocking (finally!!!) 12/2023 - Megan addd in ability to
% freeze elements 1/2/24 - MAS fixed solution guess subroutines to use T
% dependent Eg for finding Ef gueses.  Cleaned up code. 1/3/24 - MAS
% editing of fixed elements to speed up, correction of denominators in site
% blocking (1-sum not sum-1 for the finite case). 10-4-2024 added option to
% supply a vector of initial guesses for Ef_mu_vec that will be helpful for
% generalized quenching. 10-3-2025 - Aadi added output of formation
% energies G for cs and defects, Mike cleaned it up some 10-6-25: renamed
% as defect_equilibrium_dark.m (dropping the "with frozen" from the end)
% Mike simplified N_chargestates, removing duplicate blocks from the
% if/elseif options 10/10/25
% improvements 10/11/2025 by MS: 
% simplified repeated blocks of code in chargestate number calculations
% added flag to tell particleswarm to run quietly 
% updated definition of the quantum vibent per mode to (n+1)ln(n+1)-n ln(n)
% added dG for defects and chargestates as output variables for both equilibrium and fullquench
% reordered subroutines so the big long one with ten-deep for loops is at
% the bottom so you dont scroll past it any more
% standardized variable names to dummy_X (and got rid of misspelled dumy_X)


function [equilib_dark_sol] = defect_equilibrium_dark(dummy_conditions, dummy_defects)


%%%%%%% build up holder variables for outputs %%%%%%%%%%
N_chargestates_out = zeros(dummy_conditions.num_T_equilibrium,dummy_defects.num_chargestates);   % in the final output, each charge state has a column, each row is a different T. Need this different name for the holder matrix because of scope of n_defects in some of the functions
dG_cs_out = zeros(dummy_conditions.num_T_equilibrium,dummy_defects.num_chargestates);
N_defects_out = zeros(dummy_conditions.num_T_equilibrium,dummy_defects.num_defects);
dG_defects_out = zeros(dummy_conditions.num_T_equilibrium,dummy_defects.num_defects);
n_out = zeros(dummy_conditions.num_T_equilibrium,1);   %holder column vector # rows same as T
p_out = zeros(dummy_conditions.num_T_equilibrium,1);
sth1_out = zeros(dummy_conditions.num_T_equilibrium,1);
sth2_out = zeros(dummy_conditions.num_T_equilibrium,1);
EF_out = zeros(dummy_conditions.num_T_equilibrium,1);
mu_out = zeros(dummy_conditions.num_T_equilibrium,dummy_conditions.num_elements);
tot_bal_err_out = zeros(dummy_conditions.num_T_equilibrium,1);
charge_bal_err_out = zeros(dummy_conditions.num_T_equilibrium,1);
element_bal_err_out = zeros(dummy_conditions.num_T_equilibrium,1);



%% Tloop_conditions holds *scalar* values for all properties one at a time as the T loop is worked throguh
%%% make a local copy of the conditions variable to use with T dependnet
%%% properties.  At every T we assign a new value of each property to each
%%% variable as needed.  Future fix: assign T dependent dG0 values for each
%%% chargestate for example.
%% all the subroutines will be written to use a dummy name version of conditions as a dummy variable but when actually called should be called with "Tloop_conditions"
Tloop_conditions = dummy_conditions;

%%%% detect the type of calculation requested.  1=nothing frozen, 2=dfects
%%%% frozen, 3=elements frozen, 4 = defects and elments frozen.  The calc
%%%% type gets written back into the conditions variable
[Tloop_conditions] = Calc_Method(Tloop_conditions);


% check if guesses are supplied for Ef_mu_vec or not, then check if they
% are valid
if strcmp(dummy_conditions.Ef_mu_guesses_supplied_flag, 'On')
    if isempty(dummy_conditions.Ef_fixed_mu_guesses) || (max(size(dummy_conditions.Ef_fixed_mu_guesses)) ~= dummy_conditions.num_T_equilibrium) || (min(size(dummy_conditions.Ef_fixed_mu_guesses))~=dummy_conditions.num_fixed_elements + 1)
        error('If the conditions.Ef_fixed_mu_guesses field is present, it must have same number of rows as T_equilibrium and columns as num_fixed_elements + 1 and all be non-empty');
    elseif (max(size(dummy_conditions.Ef_fixed_mu_guesses)) == dummy_conditions.num_T_equilibrium) && (min(size(dummy_conditions.Ef_fixed_mu_guesses))==dummy_conditions.num_fixed_elements+1)
        disp('Ef_fixed_mu_guesses have been supplied.  It is recommended to use this only if you really need to override the automatic guess generation and temperature dependent solution following...')
    else
        error('something strange about the supplied conditions.Ef_fixed_mu_guesses')
    end
end




%%%%%%%%% Main calculation - looped over all the equilibrium T's given.
%%%%%%%%% Tloop_conditions holds only scalar variables and is updated for
%%%%%%%%% each temperature, leaving the original conditions variabl alone
disp('as of now, T dependences for mu and band parameters are being used.  Must add capabilities for other properties being T dependent if desired in future')


%% start the main loop over temperatures
for i1 = 1:Tloop_conditions.num_T_equilibrium

    Tloop_conditions.T_equilibrium = dummy_conditions.T_equilibrium(i1);  % each time through the T loop, set the values of T_equilibrium
    Tloop_conditions.kBT_equilibrium = dummy_conditions.kBT_equilibrium(i1);
    Tloop_conditions.num_T_equilibrium = 1;
    disp(strcat('Starting calc for T=',num2str(Tloop_conditions.T_equilibrium),'K'))  % call out the start of calc for each T

    % set the mu values for the current temperature.
    Tloop_conditions.muT_equilibrium = dummy_conditions.muT_equilibrium(i1,:);  % the mu values set here will be overwritten (ignored) for fixed elements as needed further down in the code

    % if Ef_fixed_mu_guesses are supplied in the conditions variable
    if strcmp(Tloop_conditions.Ef_mu_guesses_supplied_flag,'On')
        Tloop_conditions.Ef_fixed_mu_guesses = dummy_conditions.Ef_fixed_mu_guesses(i1,:);  % set the values to the i1'st row of the array
    end

    % Handle set up the ability to specify fixed defect concentrations that
    % vary vs T
    if strcmp(dummy_conditions.T_dep_fixed_defect_flag,'On')
        if size(dummy_conditions.fixed_defects_concentrations,1)==Tloop_conditions.num_T_equilibrium
            Tloop_conditions.fixed_defects_concentrations = dummy_conditions.fixed_defects_concentrations(i1,:);          % take the i1 row of the fixed defects variable from the original conditions file and assign that row in the Tloop_conditions
            Tloop_conditions.fixed_defects_concentrations = dummy_conditions.fixed_defects_concentrations';  % rotate this to a column vector (doing it in 2 steps to make it obvious we did it)
        else
            error('If T_dep_fixed_defect_flag is ON, then conditions.fixed_defect_concentrations has to be a vector with an entry for each T')
        end
    elseif strcmp(dummy_conditions.T_dep_fixed_defect_flag,'Off')
        if size(dummy_conditions.fixed_defects_concentrations,2)~=1
            error('If T_dep_fixed_defect_flag is OFF, then conditions.fixed_defect_concentrations has to be a single value (scalar) appropriate for all temperatures')
        end
    end

    % handle T dependent dH values here
    %     if strcmp(dummy_conditions.T_dep_dH,'On')
    %         Tloop_defects.dH = dummy_defects.EformT(i1,:);
    %
    %     elseif strcmp(dummy_conditions.T_dep_dH,'Off')
    %         %% do nothing different
    %     else
    %         error('T dependent fixed defects and or T dependent Eform flags
    %         must be on or off')
    %     end
    %

    % Handle T dependent Nc, Nv, and band edges
    if strcmp(dummy_conditions.T_dep_bands_flag,'On')
        Tloop_conditions.Eg = dummy_conditions.EgT_equilibrium(i1);  % the RHS of these calls need to call into original conditions variable (in case replace all is used by accident)
        Tloop_conditions.Ec = dummy_conditions.EcT_equilibrium(i1);
        Tloop_conditions.Ev = dummy_conditions.EvT_equilibrium(i1);
        Tloop_conditions.Nc = dummy_conditions.NcT_equilibrium(i1);
        Tloop_conditions.Nv = dummy_conditions.NvT_equilibrium(i1);
    elseif strcmp(dummy_conditions.T_dep_bands_flag,'Off')
        Tloop_conditions.Eg = dummy_conditions.EgRef;
        Tloop_conditions.Ec = dummy_conditions.EcRef;
        Tloop_conditions.Ev = 0;
        Tloop_conditions.Nc = dummy_conditions.NcRef;
        Tloop_conditions.Nv = dummy_conditions.NvRef;
    else
        error('T dependent band parameters must be ON or OFF')
    end



    %% run the main calculation according to the method (method won't change vs T)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% These are cases with no elements frozen. 1 is no defects frozen, 2 is some defects frozen.  We get the mu values from the conditions directly
    if Tloop_conditions.calc_method==1 || Tloop_conditions.calc_method ==2

        disp('Calculation Method 1 or 2 running please be patient....')

        if strcmp(Tloop_conditions.Ef_mu_guesses_supplied_flag,'Off')
            if i1==1
                disp('Starting grid search of Ef to get close to first solution')
                EF_guess = EF_Guess12(Tloop_conditions);   % get a guess for EF close to minimum
                disp('Found guess for first T, passing it off to fzero')
            elseif i1==2
                EF_guess = EF_from_last_T;
                disp('Using solution from prior temperature as guess for solution - make sure to space T values close enough for continuity')
            elseif i1>=3
                EF_guess = EF_from_last_T + 0.5*(Tloop_conditions.T_equilibrium - last_T)*(EF_from_last_T - EF_from_2nd_last_T)/(last_T - second_last_T);
                disp('Using avg of prior solution and 1st order Taylor expansion as guess...')
            end
        elseif strcmp(Tloop_conditions.Ef_mu_guesses_supplied_flag,'On')
            EF_guess = Tloop_conditions.Ef_mu_guuesses;
        else
            error('problem with Ef_fixed_mu_guesses - it exists but something wrong')
        end

        EF_out(i1) = fzero(@Charge_Balance12,EF_guess);
        [n_out(i1), p_out(i1), sth1_out(i1), sth2_out(i1)] = Carrier_Concentrations(EF_out(i1));  %%% carriers
        EF_full_mu_vec = [EF_out(i1) Tloop_conditions.muT_equilibrium];   % create the full EF_mu_vec taking in the fixed mu values
        [N_chargestates_out(i1,:), dG_cs_out(i1,:)] = Chargestate_Concentrations(EF_full_mu_vec);  % compute the chargestate and carrier concentrations from the full EF_mu_vec
        charge_bal_err_out(i1) = Charge_Balance12(EF_out(i1));   % compute the net charge just to check final answers
        mu_out(i1,:) = Tloop_conditions.muT_equilibrium;

        if i1>=2
            EF_from_2nd_last_T = EF_from_last_T;
            EF_from_last_T = EF_out(i1);
            second_last_T = last_T;
            last_T = Tloop_conditions.T_equilibrium;
        elseif i1==1
            EF_from_last_T = EF_out(i1);
            last_T = Tloop_conditions.T_equilibrium;
        else
            error('Somethign wrong with i1 indexing')
        end


        %% These are the cases where some elements are frozen.  3 is only elements frozen, 4 is elements and defects frozen
    elseif Tloop_conditions.calc_method==3 || Tloop_conditions.calc_method==4
        disp('Calculation Method 3 or 4 running please be patient....')
        disp('WARNING!!!: Remember, when element concentrations are fixed, it is possible that the desired solution is impossible to find within the limits of mu the user specified for each fixed element.  So check the results carefully and spend some time bracketing the solution by hand.  Remember that only complexes couple impurity elements together so you can vary one mu at a time for the most part. ')

        if strcmp(Tloop_conditions.search_method_flag,'particleswarm_pattern') || strcmp(Tloop_conditions.search_method_flag,'particleswarm_pattern_simplex')
            nvars = Tloop_conditions.num_fixed_elements + 1;
            iter_lim = Tloop_conditions.fixed_elements_swarm_iterlimit_base*nvars;
            stall_lim = Tloop_conditions.fixed_elements_swarm_stall_lim;
            display_int = Tloop_conditions.fixed_elements_swarm_display_int;

            if i1 == 1
                EF_min = Tloop_conditions.Ev - 5*Tloop_conditions.kBT_equilibrium;
                EF_max = Tloop_conditions.Ec + 5*Tloop_conditions.kBT_equilibrium;
                lb = [EF_min; Tloop_conditions.fixed_elements_mu_ranges(Tloop_conditions.indices_of_fixed_elements,1)];
                ub = [EF_max; Tloop_conditions.fixed_elements_mu_ranges(Tloop_conditions.indices_of_fixed_elements,2)];
                particles_per_fixed_mu = Tloop_conditions.fixed_elements_swarm1_fine_factor * (ub-lb) / Tloop_conditions.kBT_equilibrium;
                grid_product = prod(particles_per_fixed_mu);  % number of points needed to fill the search hypervolume spaced kBT apart
                swarmsize1 = ceil( max(Tloop_conditions.fixed_elements_swarm1_min_size, min(Tloop_conditions.fixed_elements_swarm1_max_size, grid_product^((nvars-1)/nvars)))); % initiate a very large swarm to densely sample the space.  the idea is that you want to scale better than the dimensionality of the space being searched
                if strcmp(Tloop_conditions.Ef_mu_guesses_supplied_flag,'Off')
                    if strcmp(Tloop_conditions.search_method_quiet_flag,'quiet')
                        Tloop_conditions.swarm_options1 = optimoptions('particleswarm','SwarmSize',swarmsize1,'ObjectiveLimit',Tloop_conditions.fixed_elements_swarm1_ObjectiveLimit,'MaxIterations',iter_lim,'MaxStallIterations',stall_lim,'HybridFcn', @patternsearch);
                    elseif strcmp(Tloop_conditions.search_method_quiet_flag,'verbose')
                        Tloop_conditions.swarm_options1 = optimoptions('particleswarm','SwarmSize',swarmsize1,'Display','iter','DisplayInterval',display_int,'ObjectiveLimit',Tloop_conditions.fixed_elements_swarm1_ObjectiveLimit,'MaxIterations',iter_lim,'MaxStallIterations',stall_lim,'HybridFcn', @patternsearch);
                    else
                        error('Tloop_conditions.search_method_quiet_flag must be quiet or verbose')
                    end
                elseif strcmp(Tloop_conditions.Ef_mu_guesses_supplied_flag,'On')
                    if strcmp(Tloop_conditions.search_method_quiet_flag,'quiet')
                        Tloop_conditions.swarm_options1 = optimoptions('particleswarm','InitialPoints',Tloop_conditions.Ef_fixed_mu_guesses(2:end,:),'SwarmSize',swarmsize1,'ObjectiveLimit',Tloop_conditions.fixed_elements_swarm1_ObjectiveLimit,'MaxIterations',iter_lim,'MaxStallIterations',stall_lim,'HybridFcn', @patternsearch);
                    elseif strcmp(Tloop_conditions.search_method_quiet_flag,'verbose')
                        Tloop_conditions.swarm_options1 = optimoptions('particleswarm','InitialPoints',Tloop_conditions.Ef_fixed_mu_guesses(2:end,:),'SwarmSize',swarmsize1,'Display','iter','DisplayInterval',display_int,'ObjectiveLimit',Tloop_conditions.fixed_elements_swarm1_ObjectiveLimit,'MaxIterations',iter_lim,'MaxStallIterations',stall_lim,'HybridFcn', @patternsearch);
                    else
                        error('Tloop_conditions.search_method_quiet_flag must be quiet or verbose')
                    end
                else
                    error('looks like something contradictory with Ef_fixed_mu_guesses being supplied but no valid value provided?')
                end
                EF_fixed_mu_vec_sol = particleswarm(@Total_Balance34,nvars,lb,ub,Tloop_conditions.swarm_options1);
                if strcmp(Tloop_conditions.search_method_flag,'particleswarm_pattern_simplex')
                    fminsearch_options = optimset('MaxFunEvals',Tloop_conditions.fixed_elements_fmin_MaxFunEvals,'MaxIter',Tloop_conditions.fixed_elements_fmin_MaxIter,'TolX',Tloop_conditions.fixed_elements_fmin_TolX,'TolFun',Tloop_conditions.fixed_elements_fmin_TolFun);  % may need to tune these settings
                    EF_fixed_mu_vec_sol = fminsearch(@Total_Balance34,EF_fixed_mu_vec_sol, fminsearch_options);
                end

            elseif i1>=2
                lb = EF_fixed_mu_from_last_T - Tloop_conditions.fixed_elements_swarm2_search_band_kB*Tloop_conditions.kBT_equilibrium;
                ub = EF_fixed_mu_from_last_T + Tloop_conditions.fixed_elements_swarm2_search_band_kB*Tloop_conditions.kBT_equilibrium;
                particles_per_fixed_mu = Tloop_conditions.fixed_elements_swarm2_fine_factor*Tloop_conditions.fixed_elements_swarm2_search_band_kB;
                grid_product = prod(particles_per_fixed_mu);  % number of points needed to fill the search hypervolume spaced kBT apart
                swarmsize2 = ceil(max(Tloop_conditions.fixed_elements_swarm2_min_size, min(Tloop_conditions.fixed_elements_swarm2_max_size, grid_product^((nvars-1)/nvars)))); % initiate a very large swarm to densely sample the space.  the idea is that you want to scale better than the dimensionality of the space being searched
                if strcmp(Tloop_conditions.Ef_mu_guesses_supplied_flag,'Off')
                    if strcmp(Tloop_conditions.search_method_quiet_flag,'quiet')
                        Tloop_conditions.swarm_options2 = optimoptions('particleswarm','SwarmSize',swarmsize2,'ObjectiveLimit',Tloop_conditions.fixed_elements_swarm2_ObjectiveLimit,'MaxIterations',iter_lim,'MaxStallIterations',stall_lim,'HybridFcn', @patternsearch);
                    elseif strcmp(Tloop_conditions.search_method_quiet_flag,'verbose')
                        Tloop_conditions.swarm_options2 = optimoptions('particleswarm','SwarmSize',swarmsize2,'Display','iter','DisplayInterval',display_int,'ObjectiveLimit',Tloop_conditions.fixed_elements_swarm2_ObjectiveLimit,'MaxIterations',iter_lim,'MaxStallIterations',stall_lim,'HybridFcn', @patternsearch);
                    else
                        error('Tloop_conditions.search_method_quiet_flag must be quiet or verbose')
                    end
                elseif strcmp(Tloop_conditions.Ef_mu_guesses_supplied_flag,'On')
                    if strcmp(Tloop_conditions.search_method_quiet_flag,'quiet')
                        Tloop_conditions.swarm_options2 = optimoptions('particleswarm','InitialPoints',Tloop_conditions.Ef_fixed_mu_guesses(2:end,:),'SwarmSize',swarmsize1,'ObjectiveLimit',Tloop_conditions.fixed_elements_swarm2_ObjectiveLimit,'MaxIterations',iter_lim,'MaxStallIterations',stall_lim,'HybridFcn', @patternsearch);
                    elseif strcmp(Tloop_conditions.search_method_quiet_flag,'verbose')
                        Tloop_conditions.swarm_options2 = optimoptions('particleswarm','InitialPoints',Tloop_conditions.Ef_fixed_mu_guesses(2:end,:),'SwarmSize',swarmsize1,'Display','iter','DisplayInterval',display_int,'ObjectiveLimit',Tloop_conditions.fixed_elements_swarm2_ObjectiveLimit,'MaxIterations',iter_lim,'MaxStallIterations',stall_lim,'HybridFcn', @patternsearch);
                    else
                        error('Tloop_conditions.search_method_quiet_flag must be quiet or verbose')
                    end
                else
                    error('looks like something contradictory with Ef_fixed_mu_guesses being supplied but no valid value provided?')
                end
                EF_fixed_mu_vec_sol = particleswarm(@Total_Balance34,nvars,lb,ub,Tloop_conditions.swarm_options2);
                if strcmp(Tloop_conditions.search_method_flag,'particleswarm_pattern_simplex')
                    fminsearch_options = optimset('MaxFunEvals',Tloop_conditions.fixed_elements_fmin_MaxFunEvals,'MaxIter',Tloop_conditions.fixed_elements_fmin_MaxIter,'TolX',Tloop_conditions.fixed_elements_fmin_TolX,'TolFun',Tloop_conditions.fixed_elements_fmin_TolFun);  % may need to tune these settings
                    EF_fixed_mu_vec_sol = fminsearch(@Total_Balance34,EF_fixed_mu_vec_sol, fminsearch_options);
                end
            end


        elseif strcmp(Tloop_conditions.search_method_flag,'grid_fminsearch')
            if strcmp(Tloop_conditions.Ef_mu_guesses_supplied_flag,'Off')
                max_mu_range = max(abs(Tloop_conditions.fixed_elements_mu_ranges(:,2)-Tloop_conditions.fixed_elements_mu_ranges(:,1)));
                start_kBT = dummy_conditions.kBT_equilibrium(1);  %this line calls back to the original conditions variable not Tloop one.  Otherwise the search grid would become really really fine
                fine_factor = 8;
                Tloop_conditions.fixed_elements_mu_npoints = ceil(fine_factor*max_mu_range/start_kBT);
                fminsearch_options = optimset('MaxFunEvals',Tloop_conditions.fixed_elements_fmin_MaxFunEvals,'MaxIter',Tloop_conditions.fixed_elements_fmin_MaxIter,'TolX',Tloop_conditions.fixed_elements_fmin_TolX,'TolFun',Tloop_conditions.fixed_elements_fmin_TolFun);  % may need to tune these settings
                if i1==1
                    disp('Starting grid search over Ef and mu for fixed elements to get close to the solution for the first T.  Remember that the number of grid points to look at scales as num^dim+1 where dim is the number of fixed elements and the +1 is for Ef');
                    tic
                    EF_fixed_mu_guess = EF_Fixed_Mu_Guess34(Tloop_conditions); % get a guess for EF and the fixed mu values close to minimum using trial and error
                    guess_time = toc;
                    disp(strcat('Found guess for EF and mu values for fixed elements for first T in_',num2str(guess_time),'_sec.  Phew that took a while.  I need a shower, a beer, and to put my feet up.  Wait, you need me to keep going?  Fine if you insist, but you should try to bracket the solution better next time.'));
                elseif i1==2
                    EF_fixed_mu_guess = EF_fixed_mu_from_last_T;
                    disp('Using solution from prior temperature as guess for solution - make sure to space T values close enough for continuity')
                elseif i1>=3
                    % this next lines calls back to the original conditions
                    % variable so it can see all the temperatures
                    EF_fixed_mu_guess = EF_fixed_mu_from_last_T + 0.5*((EF_fixed_mu_from_last_T - EF_fixed_mu_from_2nd_last_T)/(dummy_conditions.T_equilibrium(i1-1) - dummy_conditions.T_equilibrium(i1-2)))*(dummy_conditions.T_equilibrium(i1) - dummy_conditions.T_equilibrium(i1-1));
                    disp('Using avg of prior solution and 1st order Taylor expansion as guess...')
                end

            elseif strcmp(Tloop_conditions.Ef_mu_guesses_supplied_flag,'On')
                EF_fixed_mu_guess = Tloop_conditions.Ef_fixed_mu_guesses;
            else
                error('looks like something contradictory with Ef_fixed_mu_guesses being supplied but no valid value provided?')
            end

            % do the optimization using fminsearch.
            [EF_fixed_mu_vec_sol,~,fminsearch_exit_flag] = fminsearch(@Total_Balance34, EF_fixed_mu_guess, fminsearch_options);  % up to here we deal with EF_fixed_mu_vec.  After the optimization we rebuild the whole EF_full_mu_vec
            [EF_fixed_mu_vec_sol,~,fminsearch_exit_flag] = fminsearch(@Total_Balance34, EF_fixed_mu_vec_sol, fminsearch_options);
            % toc

        else
            error('For case 3 & 4 calculations, search_method_flag must be one of the coded options')
        end


        % these items here are executed for cases 3&4 using any solver
        % option
        disp(strcat('Found solution as__',num2str(EF_fixed_mu_vec_sol),'__for EF_fixed_mu_vec'))
        EF_full_mu_vec_sol = Expand_Fixed_Mu_Vec_To_Full_Mu_Vec(EF_fixed_mu_vec_sol);  % expand it out to the full EF_mu_vec
        EF_out(i1) = EF_full_mu_vec_sol(1);  % pick off the optimized EF value
        mu_out(i1,:) = EF_full_mu_vec_sol(2:end);
        [n_out(i1), p_out(i1), sth1_out(i1), sth2_out(i1)]  = Carrier_Concentrations(EF_out(i1));  %%% carriers
        [N_chargestates_out(i1,:),dG_cs_out(i1,:)] = Chargestate_Concentrations(EF_full_mu_vec_sol);  %% chargestates
        charge_bal_err_out(i1) = Charge_Balance34(EF_fixed_mu_vec_sol);
        element_bal_err_out(i1) = sum(Fixed_Element_Balance34(EF_fixed_mu_vec_sol));
        tot_bal_err_out(i1) = Total_Balance34(EF_fixed_mu_vec_sol);
        disp(strcat('charge_bal=',num2str(charge_bal_err_out(i1))));
        disp(strcat('element_bal=',num2str(element_bal_err_out(i1))));
        disp(strcat('tot_bal=',num2str(tot_bal_err_out(i1))));


        if i1>=2
            EF_fixed_mu_from_2nd_last_T = EF_fixed_mu_from_last_T;
            EF_fixed_mu_from_last_T = EF_fixed_mu_vec_sol;
            second_last_T = last_T;
            last_T = Tloop_conditions.T_equilibrium;
        elseif i1==1
            EF_fixed_mu_from_last_T = EF_fixed_mu_vec_sol;
            last_T = Tloop_conditions.T_equilibrium;
        else
            error('Somethign wrong with i1 indexing')
        end

    else
        error('calc method must be defined as 1,2,3, or 4 - check problem definition / setup')
    end % close the if/elseif loop for different types of calcualtions

end % close the for loop over temperatures indexed by i1

%% %%%%%%%%%%%% end of main calc %%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%% create final output strcutrues  %%%%%%%%%%%%%%%%%%%%%%%%
%% now that we are out of the T loop we go back to using dummy_conditions

%% Compute defect numbers by summing over chargstates in each defect.  Do this here at the end so the loop only runs one time %%
defect_prefac = zeros(dummy_conditions.num_T_equilibrium, dummy_defects.num_defects);
for i2 = 1:dummy_defects.num_defects
    indices_of_all_cs_in_defect = find(dummy_defects.cs_defect_ID == i2);  % will produce a column vector as long as cs_defect_ID is a column vector
    N_defects_out(:,i2) = sum(N_chargestates_out(:,indices_of_all_cs_in_defect),2);  %sum all the chargestates in the defect to get its total number (col vector)
    for i3 = 1:dummy_conditions.num_T_equilibrium
        [~,dominant_cs_index] = max(N_chargestates_out(i3, indices_of_all_cs_in_defect));  %figure out the realtive position (1st, 2nd, 3rd etc) in the defect's chargestate list of the defect with highest concentration at each T
        defect_prefac(i3, i2) = dummy_defects.cs_prefactor(indices_of_all_cs_in_defect(dominant_cs_index(1)));    % store the concentration prefactor for that dominant chargestate of the current defect, at each temperature.  %NOTE: there could be a low-probability case where the max() gives more than one index out - if two chargestates have exactly the same concentration.  what we have done is just taken the index of the first one in the list.
    end
end
% compute the effective formation energy for each defect.
dG_defects_out = -(dummy_conditions.kBT_equilibrium' * ones(1,dummy_defects.num_defects)) .* log(N_defects_out ./ defect_prefac);   % calulate the effective formation energy of the defect as a whole.  For each T and each defect, use a potentially different prefactor.  Why?  The dominant chargestate could swap vs T from one that is on a lattice site to another that has a higher geometric degeneracy factor, for example.


clear indices_of_all_cs_in_defect defect_prefac dominant_cs_index i2 i3


%% build up the solution structure
equilib_dark_sol.defect_names = dummy_defects.defect_names;
equilib_dark_sol.chargestate_names = dummy_defects.chargestate_names;
equilib_dark_sol.T_equilibrium = (dummy_conditions.T_equilibrium)';
equilib_dark_sol.Nd = dummy_conditions.Nd*ones(dummy_conditions.num_T_equilibrium,1);
equilib_dark_sol.Na = dummy_conditions.Na*ones(dummy_conditions.num_T_equilibrium,1);
equilib_dark_sol.EFn = EF_out;   % planning ahead for light calcs
equilib_dark_sol.EFp = EF_out;
equilib_dark_sol.n = n_out;
equilib_dark_sol.p = p_out;
equilib_dark_sol.sth1 = sth1_out;
equilib_dark_sol.sth2 = sth2_out;
equilib_dark_sol.chargestates = N_chargestates_out;
equilib_dark_sol.defects = N_defects_out;
equilib_dark_sol.dG_cs = dG_cs_out;
equilib_dark_sol.dG_defects = dG_defects_out;
equilib_dark_sol.charge_bal_err = charge_bal_err_out;
equilib_dark_sol.element_bal_err = element_bal_err_out;
equilib_dark_sol.tot_bal_err = tot_bal_err_out;
equilib_dark_sol.mu = mu_out;

[equilib_dark_sol.stoich, equilib_dark_sol.element_totals] = calc_stoich(equilib_dark_sol, dummy_conditions, dummy_defects);  % use the function calc_stoich.m to calculate the stoichiometry for the material of current interest.  Add new materials as files, or directly in calc_stoich.m as user pleases

%%%%%%%%% end of main function actions (close with end at end of file %%%%%%%%%%%%%%%







%% subroutines - note these all share memory space with the main function ("end" for the main function comes after all of these subroutines are ended %%


%% These functions are involved in calculating the numbers of things.  Some take in reduced and others full EF_mu_vecs holding just the fixed mu's or the full list of mu's respectively

    function [charge_bal12] = Charge_Balance12(EF_dummy)
        % function to calculate signed charge balance for methods 1 and 2
        % for speedup.  Same as the 34 version but stripped all the fixed
        % mu stuff out.  This gets called a lot so dont put any logic
        % checks inside - do checking outside before calling this
        [n, p, sth1, sth2] = Carrier_Concentrations(EF_dummy);  %deal with carriers first
        EF_full_mu_vec_dummy = [EF_dummy Tloop_conditions.muT_equilibrium];   % create the full EF_mu_vec
        [N_chargestates,~] = Chargestate_Concentrations(EF_full_mu_vec_dummy);
        charge_bal12 = sum(dummy_defects.cs_charge.* N_chargestates') + p + sth1 + sth2 - n + Tloop_conditions.Nd - Tloop_conditions.Na  ;              % absolute signed charge bal comes out signed +/-
        % charge_bal12 = (sum(dummy_defects.cs_charge.* N_chargestates') +
        % p - n + Tloop_conditions.Nd -
        % Tloop_conditions.Na)/sqrt(Tloop_conditions.Nc*Tloop_conditions.Nv)
        % ;  % scaled to ni2 sort of.   comes out signed +/-
    end %charge_bal12



    function [charge_bal34] = Charge_Balance34(EF_fixed_mu_vec_dummy)
        % function to calculate signed charge balance.  This gets called a
        % lot so dont put any logic checks inside - do checking outside
        % before calling this
        [n, p, sth1, sth2] = Carrier_Concentrations(EF_fixed_mu_vec_dummy(1));  % deal with carriers first
        % now deal with defects.  Need to take in reduced mu vec and expand
        % out to the full one before sending to chargestate calculation
        [EF_full_mu_vec_dummy] = Expand_Fixed_Mu_Vec_To_Full_Mu_Vec(EF_fixed_mu_vec_dummy);
        [N_chargestates,~] = Chargestate_Concentrations(EF_full_mu_vec_dummy);  % send this full vec to calculate the chargestates
        charge_bal34 = sum(dummy_defects.cs_charge.* N_chargestates') + p + sth1 + sth2 - n + Tloop_conditions.Nd - Tloop_conditions.Na;   % absolute signed charge bal
        % charge_bal34 = (sum(dummy_defects.cs_charge.* N_chargestates') +
        % p - n + Tloop_conditions.Nd -
        % Tloop_conditions.Na)/sqrt(Tloop_conditions.Nc*Tloop_conditions.Nv);
        % % relative comes out signed +/-
    end %end charge_bal34



    function [n, p, sth1, sth2] = Carrier_Concentrations(EF_dummy)
        %%%% function to compute n and p.  Input must be EF only (scalar
        %%%% not vector with mu's attached)  %%%%%%%%%%%%
        etaCB = (EF_dummy - Tloop_conditions.Ec)/Tloop_conditions.kBT_equilibrium;
        etaVB = (EF_dummy - Tloop_conditions.Ev)/Tloop_conditions.kBT_equilibrium;   % this looks wrong (not symmetric compared to CB case) but it is right. Direction of integration and sign on EF-Ev are both swapped.
        etaSTH1 = (EF_dummy - (Tloop_conditions.Ev + Tloop_conditions.E_relax_sth1) )/Tloop_conditions.kBT_equilibrium;
        etaSTH2 = (EF_dummy - (Tloop_conditions.Ev + Tloop_conditions.E_relax_sth2) )/Tloop_conditions.kBT_equilibrium;

        if strcmp(Tloop_conditions.Boltz_or_FD_flag,'FD')            % use Fermi-Dirac integrals so degenerate conditions handled correctly
            n = n_Fermi_Dirac(etaCB, Tloop_conditions.Nc);
            p = n_Fermi_Dirac(-etaVB, Tloop_conditions.Nv);
            sth1 = Tloop_conditions.sth_flag * n_Fermi_Dirac(-etaSTH1, Tloop_conditions.num_sites(3));  % number of STH has prefactor for all the oxygen sites that can support them
            sth2 = Tloop_conditions.sth_flag * n_Fermi_Dirac(-etaSTH2, Tloop_conditions.num_sites(4));
        elseif strcmp(Tloop_conditions.Boltz_or_FD_flag,'Boltz')              % %     use just Boltzmann approx
            n = Tloop_conditions.Nc * exp(etaCB);   % these are right (sign swap).  Boltzmann factors should end up <1 when EF is in gap
            p = Tloop_conditions.Nv * exp(-etaVB);
            sth1 = Tloop_conditions.sth_flag * Tloop_conditions.num_sites(3) * exp(-etaSTH1);
            sth2 = Tloop_conditions.sth_flag * Tloop_conditions.num_sites(4) * exp(-etaSTH2);
        else
            error('Boltz_or_FD_flag must be Boltz or FD')
        end
    end



    function [tot_bal34] = Total_Balance34(EF_fixed_mu_vec_dummy)
        %% this constructs the objective function that is to be zeroed out.  Takes in the reduced EF_mu_Vec and gives a scalar out
        % This is the log10 of the total error from zero on the element sum
        % and charge sum
        %         tot_bal34 =
        %         log10(abs(Charge_Balance34(EF_fixed_mu_vec_dummy)) +
        %         sum(abs(Fixed_Element_Balance34(EF_fixed_mu_vec_dummy)),'omitnan'));
        %         tot_bal34 =
        %         log10(abs(Charge_Balance34(EF_fixed_mu_vec_dummy)) +
        %         sum(abs(Fixed_Element_Balance34(EF_fixed_mu_vec_dummy))));
        tot_bal34 = abs(Charge_Balance34(EF_fixed_mu_vec_dummy)) + sum(abs(Fixed_Element_Balance34(EF_fixed_mu_vec_dummy)));
    end  %%% end tot_balance objective function  %%%%%%%%%%%%%



    function [fixed_element_bal34_vec] = Fixed_Element_Balance34(EF_fixed_mu_vec_dummy)
        % returns a vector of the signed difference between fixed element
        % concs and their targets, only for the fixed elements
        %         element_bal34 =
        %         zeros(size(Tloop_conditions.num_fixed_elements));  %
        %         initialize as all zeros
        [EF_full_mu_vec_dummy] = Expand_Fixed_Mu_Vec_To_Full_Mu_Vec(EF_fixed_mu_vec_dummy);
        [N_chargestates,~] = Chargestate_Concentrations(EF_full_mu_vec_dummy);  % send this full vec to calculate the chargestates
        %         fixed_element_bal34_vec =
        %         ((((dummy_defects.cs_dm')*N_chargestates')' -
        %         Tloop_conditions.fixed_conc_values).*Tloop_conditions.fixed_conc_flag)./Tloop_conditions.fixed_conc_values;
        %         % have to rotate and reshape the arrays a few times to
        %         get sizes and shapes right.  Then .* by the fixed element
        %         mask to isolate only the target elements - for example a
        %         substitutional defect also produces a missing host atom
        %         so dont want to include those. This masking is applied
        fixed_element_bal34_vec =  (((dummy_defects.cs_dm')*N_chargestates')' - Tloop_conditions.fixed_conc_values).*Tloop_conditions.fixed_conc_flag;   % have to rotate and reshape the arrays a few times to get sizes and shapes right.  Then .* by the fixed element mask to isolate only the target elements - for example a substitutional defect also produces a missing host atom so dont want to include those. This masking is applied
    end % end element bal34



    function [N_chargestates, dG_cs_out] = Chargestate_Concentrations(EF_full_mu_vec_dummy)
        % function that calculates the numbers of each charge state of
        % each defect.  Input has to be a full EF_mu_vector with
        % entries for ALL elements not just fixed ones.  Fixed defects
        % overrides fixed elements, thus the 4 method cases are spit
        % according to whether or not any defects are fixed.  The main
        % calc loop takes care of supplying a full EF_mu_vec to send in.

        % figure out what to do about -TdS_vib term
        if strcmp(Tloop_conditions.vib_ent_flag,'3kB')
            dG_vib_ent = 3 * Tloop_conditions.kBT_equilibrium * sum(dummy_defects.cs_dm,2);  % this is TdS.  A vacacny has sum(cs_dm)=-1, and interstitial has +1.  dSvib = 3*kB*sum(cs_dm)*f(T) where f(T) is the classical or quantum function (positive numbers).  dG=dH-TdS = dH - 3kBT*f(T)*sum(dm).
        elseif strcmp(Tloop_conditions.vib_ent_flag,'Quantum')
            dG_vib_ent = 3 * Tloop_conditions.kBT_equilibrium * dSvib_quantum_per_mode(Tloop_conditions.T_equilibrium, Tloop_conditions.vibent_T0) * sum(dummy_defects.cs_dm,2);
        elseif strcmp(Tloop_conditions.vib_ent_flag,'Off')
            dG_vib_ent = 0;
        else
            error('Vibrational entropy flag must be 3kB, Quantum, or Off')
        end

        % compute total and relative formation energies (total has chemical
        %                                     
        % potential part, relative doesnt).
        % fixed dH0, then q*Ef term, then -TDSvib term.
        dG_rel = dummy_defects.cs_Eform + dummy_defects.cs_charge*EF_full_mu_vec_dummy(1) - dG_vib_ent;
        dG_cs_out = dG_rel - dummy_defects.cs_dm*EF_full_mu_vec_dummy(2:end)';  % For the full dG the last thing is mu times number of atoms

        % compute concentrations from dG depending on what kind of site statistics chosen
        % Note: we compute them all here, and adjust the ones we need below
        %
        if strcmp(Tloop_conditions.site_blocking_flag,'On_Infinite')
            N_chargestates = (dummy_defects.cs_prefactor .* exp(-dG_cs_out/Tloop_conditions.kBT_equilibrium) ./ (1+sum(exp(-dG_cs_out/Tloop_conditions.kBT_equilibrium)) ))';    %% with site blocking for infinite lattice
        elseif strcmp(Tloop_conditions.site_blocking_flag,'On_Finite')
            N_chargestates = (dummy_defects.cs_prefactor .* exp(-dG_cs_out/Tloop_conditions.kBT_equilibrium) ./ (1-sum(exp(-dG_cs_out/Tloop_conditions.kBT_equilibrium)) ))';    %% with site blocking for finite crystal
        elseif strcmp(Tloop_conditions.site_blocking_flag,'Off')
            N_chargestates = (dummy_defects.cs_prefactor .* exp(-dG_cs_out/Tloop_conditions.kBT_equilibrium))';   % Boltzmann dilute defects approximation here for all charge states of all defects.
        else
            error('Site blocking flag must be On_Infinite, On_Finite, or Off')
        end

        % for calc methods 1 and 3, no defects are frozen so we have what we need already so do nothing more
        % 1 is no defects or elements frozen, 3 is elments but not defects frozen.
        if Tloop_conditions.calc_method==1 || Tloop_conditions.calc_method==3

            % do nothing more

            % % These are cases with elements frozen. 2 = defects, but not elements frozen, 4 = defects and elements frozen.
        elseif Tloop_conditions.calc_method==2 || Tloop_conditions.calc_method==4
            Boltz_fac_rel = exp(-dG_rel/Tloop_conditions.kBT_equilibrium); % Boltzmann factors for all chargstates - to be used in Gibbs distributions only for fixed/frozen defects
            % Here we overwrite concentrations for fixed defects and create
            % and overwrite effective dG values for them too.
            for i17 = 1:Tloop_conditions.num_fixed_defects     % loop through the fixed defects one at a time and assign their chargestates first
                which_frozen_defects_index = Tloop_conditions.fixed_defects_index(i17);
                which_frozen_chargestates_index = dummy_defects.defect_cs_lo(which_frozen_defects_index):dummy_defects.defect_cs_hi(which_frozen_defects_index);
                Z = sum(Boltz_fac_rel(which_frozen_chargestates_index));  % this is the Z for this particular i17'th defect
                N_chargestates(which_frozen_chargestates_index) = Tloop_conditions.fixed_defects_concentrations(which_frozen_defects_index)*Boltz_fac_rel(which_frozen_chargestates_index)/Z ;     % denom is a scalar.  Compute the conc of each charge state in defect i6 is Boltz factor/Z(i)*total conc  Totalconc is a scalar also - same for all chargestates of the give defetc
                dG_cs_out(which_frozen_chargestates_index) = -Tloop_conditions.kBT_equilibrium*log(N_chargestates(which_frozen_chargestates_index)'./dummy_defects.cs_prefactor(which_frozen_chargestates_index));   % compute the dG that would give the concentration we just calculated and repalce the one calculated the normal way
            end

        else
            error('must be one of four calculation methods based on fixed/open elements and defects')
        end
    end   % end chargestate_concentrations



    function [EF_full_mu_vec_out] = Expand_Fixed_Mu_Vec_To_Full_Mu_Vec(EF_fixed_mu_vec_dummy)
        % does what the name says - takes in a reduced EF_fixed_mu_vec and
        % spits out the whole one with mu values for all elements.
        fixed_mu_vec = EF_fixed_mu_vec_dummy(2:end);  %strip off EF from front so only have the fixed mu's
        full_mu_vec = Tloop_conditions.muT_equilibrium;   % initialize as those currently in the Tloop_conditions.mu - this assigns all the NON-fixed ones as well as the fixed ones
        counter = 0;
        for i18 = Tloop_conditions.indices_of_fixed_elements  % loop calling only the fixed elements
            counter = counter + 1;
            full_mu_vec(i18) = fixed_mu_vec(counter);  % write fixed mu's into the full mu vec (overwriting anything already there).
        end
        EF_full_mu_vec_out = [EF_fixed_mu_vec_dummy(1) full_mu_vec];
    end



    function [dSvib_Q_per_mode_norm] = dSvib_quantum_per_mode(T, T0)
        % T0/T is x = hbar*omega0/kBT.  So compute the T0 for the mean mode
        % from the Debeye temperature.
        % make sure to use T and not kBT!! 
        n_ph = 1/(exp(T0/T)-1);
        dSvib_Q_per_mode_norm = (1+n_ph)*log(1+n_ph) - n_ph*log(n_ph) +log(2); 
        % dSvib_Q_per_mode_norm = T0/(2*TempK)*coth(T0/(2*TempK)) + log(csch(T0/(2*TempK)));
    end



% function [dSvib] = Quantum_dSvib(TempK, w0_added, w0_removed)
%     % % computes the dSvib for adding some phonon modes and removing %
%     others while forming a chargestate. % w0 values should be in s-1
%     units (not cm-1). % TempK is a scalar temperature. % w0_added and
%     w0_removed are matrices with each row listing the phonon frequencies
%     % added by creating the defect.  If none, just pad with zeros.  The %
%     number of columns will be determined by the chargestate that %
%     adds/subtracts the most modes.  Modes can be repeated if it
%     adds/subtracts multiple degenerate modes - jsut list % them multiple
%     times for that row.  The number of columns need not % be equal in the
%     added and removed matrices. % % Modes are treated as dirac deltas in
%     the phDOS.  Reference is the perfect unit cell or supercell
%     appropriate for the number of cells available in the crystal for that
%     defect or complex. % Changing the number of atoms in that volume can
%     be achieved by % uneuqal numbers of modes added and subtracted.  Each
%     atom % removed should subtract 3 modes, each added should add 3.  A %
%     substitutional should have equal numbers of added and % subtracted,
%     while vacancies and interstitials and complexes may % have unequal
%     numbers. % dSvib finally output should be a column vector with rows
%     euqal to
%     % # chargestates.
%
%     hbar = 1.05457182e-34; kB = 1.380649e-23;
%
%     dSvib = zeros(dummy_defects.num_chargestates,1);
%
%     for i = 1:dummy_defects.num_chargestates
%         % for the ith defect, pick its row from the w0 matrix and take %
%         only the nonzero elements of that row.  zeros are for padding %
%         to make a proper matrix.  If the x_added_vec (or the removed %
%         one) is empty, meaning that all the entries were zero, it's not a
%         problem as the dSvib will return 0. % the nonzeros() function
%         gives a column vector so need to take % its trasnpose.
%         x_added_vec = hbar*nonzeros(w0_added(i,:))' / (kB*TempK);
%         x_removed_vec = hbar*nonzeros(w0_removed(i,:))' / (kB*TempK);
%         dSvib(i) = Dirac_Delta_w0_dSvib(x_added_vec) -
%         Dirac_Delta_w0_dSvib(x_removed_vec);
%     end
% end
%
%
%
% function [dSvib] = Dirac_Delta_w0_dSvib(x_vec)
%     % this can't handle zeros in the x_vec since the Bose distribution %
%     diverges at 0. dSvib_vec = (exp(x_vec)./(exp(x_vec)-1)) .* log(
%     (1+coth(x_vec/2))/2) - (1./(exp(x_vec)-1)) .* log(1./(exp(x_vec)-1));
%     dSvib = sum(dSvib_vec);
% end




%% These functions have the conditions variable as input, but in each one we make it a dummy name.  When actually calling these functions we want to call it with "Tloop_conditions" as input.  They are written with a dummy conditions variable to differentiate the context

    function [CM_conditions] = Calc_Method(CM_conditions)
        %% check for valid imputs for fixed and open defects and elements
        if numel(CM_conditions.fixed_defects)~=dummy_defects.num_defects  % check if each defect has a 1 or 0
            error('Tloop_conditions.fixed_defects has to be same size as num_defects')
        elseif numel(CM_conditions.fixed_conc_flag)~=dummy_defects.numelements  % check that each element part of the calculation is listed
            error('Tloop_conditions.fixed_conc_flag has to be same size as # elements')
        elseif sum(CM_conditions.fixed_defects) ~= CM_conditions.num_fixed_defects
            error('Something wrong with number of fixed defects')

        elseif sum(CM_conditions.fixed_conc_flag) ~= CM_conditions.num_fixed_elements
            error('Something wrong with number of fixed elements')
        end

        %% these logic checks figure out how calc will be done
        if sum(CM_conditions.fixed_defects)==0   % no defects fixed, all open
            CM_conditions.some_defects_fixed_flag = 0;
            CM_conditions.fixed_defects_index = [];
            CM_conditions.num_fixed_defects = 0;
            % CM_conditions.open_defects_index =
            % ones(size(CM_conditions.fixed_defects));
            % CM_conditions.num_open_defects =
            % max(size(CM_conditions.open_defects_index));
        elseif sum(CM_conditions.fixed_defects)~=0  %% some defects have fixed concentration
            CM_conditions.some_defects_fixed_flag = 1;
            CM_conditions.fixed_defects_index = find(CM_conditions.fixed_defects==1);
            CM_conditions.num_fixed_defects = max(size(CM_conditions.fixed_defects_index));
            % CM_conditions.open_defects_index =
            % find(CM_conditions.fixed_defects==0);
            % CM_conditions.num_open_defects =
            % max(size(CM_conditions.open_defects_index));
        end

        if sum(CM_conditions.fixed_conc_flag)==0   % no elements fixed, all open (mu set for all not concentration)
            CM_conditions.some_elements_fixed_flag = 0;
            CM_conditions.fixed_elements_index = [];
            CM_conditions.num_fixed_elements = 0;
            CM_conditions.open_elements_index = ones(size(CM_conditions.fixed_conc_flag));
        elseif sum(CM_conditions.fixed_conc_flag)~=0  %% at least one element concentration is set
            CM_conditions.some_elements_fixed_flag = 1;
            CM_conditions.fixed_elements_index = find(CM_conditions.fixed_conc_flag==1);
            CM_conditions.num_fixed_elements = numel(CM_conditions.fixed_elements_index);
            CM_conditions.open_elements_index = find(CM_conditions.fixed_conc_flag==0);
        end

        if CM_conditions.some_defects_fixed_flag==0 && CM_conditions.some_elements_fixed_flag==0
            CM_conditions.calc_method = 1;
        elseif CM_conditions.some_defects_fixed_flag==1 && CM_conditions.some_elements_fixed_flag==0
            CM_conditions.calc_method = 2;
        elseif CM_conditions.some_defects_fixed_flag==0 && CM_conditions.some_elements_fixed_flag==1
            CM_conditions.calc_method = 3;
        elseif CM_conditions.some_defects_fixed_flag==1 && CM_conditions.some_elements_fixed_flag==1
            CM_conditions.calc_method = 4;
        else
            error('defects and elements must be either fixed or open')
        end

    end



    function [EF_guess12] = EF_Guess12(EFG12_conditions)
        %%%% this will only be called for mthods 1 & 2. for methods 1&2,
        %%%% function that uses grid search to get an EF guess close to the
        %%%% charge balance solution.  It should return 1 value (old method
        %%%% tried to give 2 but it was overcomplicated - leave it as
        %%%% comments for now in case it is needed sometime)
        % keep kBT and Eg as local variables in this function by using
        % distinct aliases for the input arguments

        [EF_grid] = EF_Grid_Maker(EFG12_conditions);

        nn = numel(EF_grid);
        errs = zeros(1,nn);
        for i4=1:nn
            errs(i4) = Charge_Balance12(EF_grid(i4));   % methods 1&2 will only involve Ef
        end

        %%% new method for finding EF guess - just find point that is
        %%% nearby the zero and hand it off to fzero
        [~,min_index] = min(errs.^2);    % find the index of the point closest to the zero crossing
        if isscalar(min_index)   %normal case where one point is closest to the zero crossing
            EF_guess12 = EF_grid(min_index);
        elseif numel(min_index)==2  % case where two are equally close, just pick one
            EF_guess12 = EF_grid(min_index(1));
        elseif isempty(min_index) || numel(min_index)>=3  %if something is really weird plot chargebal equation vs EF and give an error
            figure(1)
            clf
            hold on
            plot(EF_grid,log10(charge_element_balance),'b-')
            errs(min_index)
            error('something strange about charge balance error vs EF - solutions may not be valid.  See figure 1')
        end
        %% old method trying to find two points for x0 that are both on the same side of the zero so it points fzero downhill - probably was overkill
        %         edge_index = find(diff(sign(errs))~=0);   %  this finds
        %         the rising/falling edge where error changes sign.  The
        %         usual case is that EF small means charge bal will be +
        %         becasue of all the holes, and EF too high will be neg.
        %         This could change if sum(size(edge_index)==[1 1])==2   %
        %         this finds the case where only one index is found where
        %         sign changes
        %             min_index = [edge_index edge_index+1];   % so we are
        %             finding the two guesses that bracket the solution one
        %             + and one - %     elseif sum(size(edge_index)~=[1
        %             1])==2 %         disp('waring: charge balance error
        %             may not be monotonic')
        %         elseif sum(size(edge_index)==[1 2])==2     % this is the
        %         normal casee where two indices are found.  The derivative
        %         points towards the zero if they ar both on the same side
        %         of zero.
        %             min_index = edge_index;
        %         else
        %             figure(1) clf hold on
        %             plot(EF_grid,log10(charge_bal),'b-') errs(min_index)
        % %             error('something strange about charge balance error
        % vs EF - solutions may not be valid.  See figure 1') %         end
        % % %         guess = EF_grid(min_index);
        %
        %         if size(guess,2)~=2
        %             error('EF guess failed to find two guesses; check for
        %             unique situation')
        %         end

    end     %%%% end EF_guess12  %%%%



    function [fixed_mu_grids] = Fixed_Mu_Ggrid_Maker(FMGM_conditions)
        % this gets called only for methods 3 &4, and gets called before
        % anything else function to create an array holding vectors of mu
        % values for each frozen element (and no others).  The fixed
        % elements will be in order in the array, but not the same absolute
        % element numbers. So if there are 6 elements, and 1, 3, and 5 are
        % frozen, this will have 3 cells, 1=1st, 3=2nd , 5=3rd
        %% future possible improvement: can just make holder a matrix directly rather than making a cell array here and turning it back into a vector when used in guess34
        %% also could move this task into the script instead and include this object in the conditions variable sent.  In that case it would be better to have it as a cell array as then some entries can be nulls rathet than zeros required if this is a matrix

        fixed_mu_grids = cell(1,FMGM_conditions.num_fixed_elements);  % create empty cell array with enough entries to hold a vector of mu values for each fixed element
        counter = 1;
        for i16 = FMGM_conditions.indices_of_fixed_elements  % loop over the i16 values indexing the frozen elements  This way only executes the loop for the frozen ones, skipping the non-frozen ones so gets rid of need to insert if elseif inside the for loop.  This will only be called in methods 3 and 4 so no need to worry about methods 1&2 when none are frozen
            mu_inc = (FMGM_conditions.fixed_elements_mu_ranges(i16,2)-FMGM_conditions.fixed_elements_mu_ranges(i16,1))/FMGM_conditions.fixed_elements_mu_npoints; % decide up front that you will divide the rnage into a fixed number of increments
            fixed_mu_grids{counter} = FMGM_conditions.fixed_elements_mu_ranges(i16,1) : mu_inc : FMGM_conditions.fixed_elements_mu_ranges(i16,2);
            counter = counter + 1;
        end
    end



    function [EF_grid] = EF_Grid_Maker(EFGM_conditions)
        % Set up a grid of EF values the grid.  Using interval 1/2 of kBT
        % guarantees you can't miss the solution whcih should be thus
        % within kB/2 of the guess
        EF_inc = EFGM_conditions.kBT_equilibrium/2;
        EF_grid = EFGM_conditions.Ev-10*EF_inc : EF_inc : ceil((EFGM_conditions.Ec-EFGM_conditions.Ev)/EF_inc)*EF_inc + 10*EF_inc;  % this makes a grid to check over the range -5kBT to Eg+5kBT, the ceil() command takes care of casee where Ec-Ev/interval is not an integer
    end



    function [EF_fixed_mu_guess34] = EF_Fixed_Mu_Guess34(EFFMG34_conditions)
        % take in conditions variable and spit out a EF_fixed_mu_vec as a
        % guess.  Grid search the m dimensional space where m = number of
        % frozen elements.  Vctor length is m+1 for EF.  This guess is sent
        % to optimizer to get mu values for the fixed ones.  As of now it
        % handles up to 10 elements being fixed - to add more jsut add more
        % cases.
        %% future improvement: get rid of the explicit sets of for loops of varying recursion.  Replace with an array with elements indexed sequentially rather than having multidimensional indexes.

        %% IMPORTANT Future thing to try: monte carlo sampling for generating the guesses - especially when multiple elements are fixed this will pay off.

        %% NOTE: it may be tempting to not search as fine a grid (default is to use the lesser of range/(kBT/2) or 200 intervals).  However the solutions especially at T lower than 1000 can be in a razor edge trough surrounded by a flat landsape on either side.  So it is worth it to do a good first search before proceeding - otherwis you get garbage.

        % handle Ef first %% Set up the grid for Ef in same way as for
        % guess12
        [EF_grid] = EF_Grid_Maker(EFFMG34_conditions);

        % create search grid for each frozen element
        [mu_grid_holder] = Fixed_Mu_Ggrid_Maker(EFFMG34_conditions);   % this creates a cell array with number of cells equal to the number of fixed elements.  In those cells are vectors of mu grids for each fixed element.

        % start going through the cases for diffrent numbers of fixed
        % elements
        tot_bal = Inf;  %initialize tot_bal as infinite (< operator will pick any finite value as less than on first iteration)

        % converting this to a switch structrue so it only evaluates one
        % case rather than going through all the elseifs, since this gets
        % called a lot


        switch EFFMG34_conditions.num_fixed_elements
            case 1

                %         if EFFMG34_conditions.num_fixed_elements == 1
                mu1_grid = cell2mat(mu_grid_holder(1));   % extract the first mu vec to search over and make it a vector
                %                 diagnstics charge_holder =
                %                 zeros(numel(EF_grid),numel(mu1_grid));
                %                 element_holder =
                %                 zeros(numel(EF_grid),numel(mu1_grid));
                %                 tot_holder =
                %                 zeros(numel(EF_grid),numel(mu1_grid));

                num_grid_points = numel(EF_grid)*numel(mu1_grid);
                disp(strcat('Starting grid search over_',num2str(num_grid_points),'_points in Ef-mu space'))
                for i5 = 1:numel(EF_grid)
                    for i6 = 1:numel(mu1_grid)
                        current_EF_fixed_mu_vec = [EF_grid(i5) mu1_grid(i6)];         % create the current EF_fixed_mu_vec
                        current_tot_bal = Total_Balance34(current_EF_fixed_mu_vec); % throw it into the total balance function
                        %                         charge_holder(i5,i6) =
                        %                         Charge_Balance34(current_EF_fixed_mu_vec);
                        %                         element_holder(i5,i6) =
                        %                         sum(Fixed_Element_Balance34(current_EF_fixed_mu_vec));
                        %                         tot_holder(i5,i6) =
                        %                         abs(charge_holder(i5,i6))
                        %                         +
                        %                         abs(element_holder(i5,i6));
                        if current_tot_bal < tot_bal    % see if it is better than prior values.  If yes, set guess to this
                            EF_fixed_mu_guess34 = current_EF_fixed_mu_vec;
                            tot_bal = current_tot_bal;
                        end
                    end
                end

            case 2
                %         elseif EFFMG34_conditions.num_fixed_elements == 2
                mu1_grid = cell2mat(mu_grid_holder(1));
                mu2_grid = cell2mat(mu_grid_holder(2));
                num_grid_points = numel(EF_grid)*numel(mu1_grid)*numel(mu2_grid);
                %                 update_progress_num =
                %                 floor(num_grid_points/10);
                disp(strcat('Starting grid search over_',num2str(num_grid_points),'_points in Ef-mu space'))
                for i5 = 1:numel(EF_grid)
                    percent_done = i5/numel(EF_grid)*100;
                    disp(strcat(num2str(percent_done),'% throguh initial grid search.  Appreciate your patience'))
                    for i6 = 1:numel(mu1_grid)
                        for i7 = 1:numel(mu2_grid)
                            current_EF_fixed_mu_vec = [EF_grid(i5) mu1_grid(i6) mu2_grid(i7)];    % create the current EF_fixed_mu_vec
                            current_tot_bal = Total_Balance34(current_EF_fixed_mu_vec);
                            if current_tot_bal < tot_bal
                                EF_fixed_mu_guess34 = current_EF_fixed_mu_vec;
                                tot_bal = current_tot_bal;
                            end
                        end
                    end
                end

            case 3
                %         elseif EFFMG34_conditions.num_fixed_elements == 3
                mu1_grid = cell2mat(mu_grid_holder(1));
                mu2_grid = cell2mat(mu_grid_holder(2));
                mu3_grid = cell2mat(mu_grid_holder(3));
                num_grid_points = numel(EF_grid)*numel(mu1_grid)*numel(mu2_grid)*numel(mu3_grid);
                disp(strcat('Starting grid search over_',num2str(num_grid_points),'_points in Ef-mu space'))
                for i5 = 1:numel(EF_grid)
                    percent_done = i5/numel(EF_grid)*100;
                    disp(strcat(num2str(percent_done),'% throguh initial grid search.  Appreciate your patience'))
                    for i6 = 1:numel(mu1_grid)
                        for i7 = 1:numel(mu2_grid)
                            for i8 = 1:numel(mu3_grid)
                                current_EF_fixed_mu_vec = [EF_grid(i5) mu1_grid(i6) mu2_grid(i7) mu3_grid(i8)];    % create the current EF_fixed_mu_vec
                                current_tot_bal = Total_Balance34(current_EF_fixed_mu_vec);
                                if current_tot_bal < tot_bal
                                    EF_fixed_mu_guess34 = current_EF_fixed_mu_vec;
                                    tot_bal = current_tot_bal;
                                end
                            end
                        end
                    end
                end

            case 4
                %         elseif EFFMG34_conditions.num_fixed_elements == 4
                mu1_grid = cell2mat(mu_grid_holder(1));
                mu2_grid = cell2mat(mu_grid_holder(2));
                mu3_grid = cell2mat(mu_grid_holder(3));
                mu4_grid = cell2mat(mu_grid_holder(4));
                num_grid_points = numel(EF_grid)*numel(mu1_grid)*numel(mu2_grid)*numel(mu3_grid)*numel(mu4_grid);
                disp(strcat('Starting grid search over_',num2str(num_grid_points),'_points in Ef-mu space'))
                for i5 = 1:numel(EF_grid)
                    percent_done = i5/numel(EF_grid)*100;
                    disp(strcat(num2str(percent_done),'% throguh initial grid search.  Appreciate your patience'))
                    for i6 = 1:numel(mu1_grid)
                        for i7 = 1:numel(mu2_grid)
                            for i8 = 1:numel(mu3_grid)
                                for i9 = 1:numel(mu4_grid)
                                    current_EF_fixed_mu_vec = [EF_grid(i5) mu1_grid(i6) mu2_grid(i7) mu3_grid(i8) mu4_grid(i9)];    % create the current EF_fixed_mu_vec
                                    current_tot_bal = Total_Balance34(current_EF_fixed_mu_vec);
                                    if current_tot_bal < tot_bal
                                        EF_fixed_mu_guess34 = current_EF_fixed_mu_vec;
                                        tot_bal = current_tot_bal;
                                    end
                                end
                            end
                        end
                    end
                end

            case 5
                %         elseif EFFMG34_conditions.num_fixed_elements == 5
                mu1_grid = cell2mat(mu_grid_holder(1));
                mu2_grid = cell2mat(mu_grid_holder(2));
                mu3_grid = cell2mat(mu_grid_holder(3));
                mu4_grid = cell2mat(mu_grid_holder(4));
                mu5_grid = cell2mat(mu_grid_holder(5));
                num_grid_points = numel(EF_grid)*numel(mu1_grid)*numel(mu2_grid)*numel(mu3_grid)*numel(mu4_grid)*numel(mu5_grid);
                disp(strcat('Starting grid search over_',num2str(num_grid_points),'_points in Ef-mu space'))
                for i5 = 1:numel(EF_grid)
                    percent_done = i5/numel(EF_grid)*100;
                    disp(strcat(num2str(percent_done),'% throguh initial grid search.  Appreciate your patience'))
                    for i6 = 1:numel(mu1_grid)
                        for i7 = 1:numel(mu2_grid)
                            for i8 = 1:numel(mu3_grid)
                                for i9 = 1:numel(mu4_grid)
                                    for i10 = 1:numel(mu5_grid)  % dont use n or p as these are carrier concentrations.  q is charge
                                        current_EF_fixed_mu_vec = [EF_grid(i5) mu1_grid(i6) mu2_grid(i7) mu3_grid(i8) mu4_grid(i9) mu5_grid(i10)];    % create the current EF_fixed_mu_vec
                                        current_tot_bal = Total_Balance34(current_EF_fixed_mu_vec);
                                        if current_tot_bal < tot_bal
                                            EF_fixed_mu_guess34 = current_EF_fixed_mu_vec;
                                            tot_bal = current_tot_bal;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end

            case 6
                %         elseif EFFMG34_conditions.num_fixed_elements == 6
                mu1_grid = cell2mat(mu_grid_holder(1));
                mu2_grid = cell2mat(mu_grid_holder(2));
                mu3_grid = cell2mat(mu_grid_holder(3));
                mu4_grid = cell2mat(mu_grid_holder(4));
                mu5_grid = cell2mat(mu_grid_holder(5));
                mu6_grid = cell2mat(mu_grid_holder(6));
                num_grid_points = numel(EF_grid)*numel(mu1_grid)*numel(mu2_grid)*numel(mu3_grid)*numel(mu4_grid)*numel(mu5_grid)*numel(mu6_grid);
                disp(strcat('Starting grid search over_',num2str(num_grid_points),'_points in Ef-mu space'))
                for i5 = 1:numel(EF_grid)
                    percent_done = i5/numel(EF_grid)*100;
                    disp(strcat(num2str(percent_done),'% throguh initial grid search.  Appreciate your patience'))
                    for i6 = 1:numel(mu1_grid)
                        for i7 = 1:numel(mu2_grid)
                            for i8 = 1:numel(mu3_grid)
                                for i9=1:numel(mu4_grid)
                                    for i10=1:numel(mu5_grid)   % dont use n or p as these are carrier concentrations.  q is charge
                                        for i11=1:numel(mu6_grid)   % dont use n or p as these are carrier concentrations.  q is charge
                                            current_EF_fixed_mu_vec = [EF_grid(i5) mu1_grid(i6) mu2_grid(i7) mu3_grid(i8) mu4_grid(i9) mu5_grid(i10) mu6_grid(i11)];    % create the current EF_fixed_mu_vec
                                            current_tot_bal = Total_Balance34(current_EF_fixed_mu_vec);
                                            if current_tot_bal < tot_bal
                                                EF_fixed_mu_guess34 = current_EF_fixed_mu_vec;
                                                tot_bal = current_tot_bal;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end

            case 7
                %         elseif EFFMG34_conditions.num_fixed_elements == 7
                mu1_grid = cell2mat(mu_grid_holder(1));
                mu2_grid = cell2mat(mu_grid_holder(2));
                mu3_grid = cell2mat(mu_grid_holder(3));
                mu4_grid = cell2mat(mu_grid_holder(4));
                mu5_grid = cell2mat(mu_grid_holder(5));
                mu6_grid = cell2mat(mu_grid_holder(6));
                mu7_grid = cell2mat(mu_grid_holder(7));
                num_grid_points = numel(EF_grid)*numel(mu1_grid)*numel(mu2_grid)*numel(mu3_grid)*numel(mu4_grid)*numel(mu5_grid)*numel(mu6_grid)*numel(mu7_grid);
                disp(strcat('Starting grid search over_',num2str(num_grid_points),'_points in Ef-mu space'))
                for i5 = 1:numel(EF_grid)
                    percent_done = i5/numel(EF_grid)*100;
                    disp(strcat(num2str(percent_done),'% throguh initial grid search.  Appreciate your patience'))
                    for i6 = 1:numel(mu1_grid)
                        for i7 = 1:numel(mu2_grid)
                            for i8 = 1:numel(mu3_grid)
                                for i9=1:numel(mu4_grid)
                                    for i10=1:numel(mu5_grid)   % dont use n or p as these are carrier concentrations.  q is charge
                                        for i11=1:numel(mu6_grid)   % dont use n or p as these are carrier concentrations.  q is charge
                                            for i12=1:numel(mu7_grid)
                                                current_EF_fixed_mu_vec = [EF_grid(i5) mu1_grid(i6) mu2_grid(i7) mu3_grid(i8) mu4_grid(i9) mu5_grid(i10) mu6_grid(i11) mu7_grid(i12)];    % create the current EF_fixed_mu_vec
                                                current_tot_bal = Total_Balance34(current_EF_fixed_mu_vec);
                                                if current_tot_bal < tot_bal
                                                    EF_fixed_mu_guess34 = current_EF_fixed_mu_vec;
                                                    tot_bal = current_tot_bal;
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end

            case 8
                %         elseif EFFMG34_conditions.num_fixed_elements == 8
                mu1_grid = cell2mat(mu_grid_holder(1));
                mu2_grid = cell2mat(mu_grid_holder(2));
                mu3_grid = cell2mat(mu_grid_holder(3));
                mu4_grid = cell2mat(mu_grid_holder(4));
                mu5_grid = cell2mat(mu_grid_holder(5));
                mu6_grid = cell2mat(mu_grid_holder(6));
                mu7_grid = cell2mat(mu_grid_holder(7));
                mu8_grid = cell2mat(mu_grid_holder(8));
                num_grid_points = numel(EF_grid)*numel(mu1_grid)*numel(mu2_grid)*numel(mu3_grid)*numel(mu4_grid)*numel(mu5_grid)*numel(mu6_grid)*numel(mu7_grid)*numel(mu8_grid);
                disp(strcat('Starting grid search over_',num2str(num_grid_points),'_points in Ef-mu space'))
                for i5 = 1:numel(EF_grid)
                    percent_done = i5/numel(EF_grid)*100;
                    disp(strcat(num2str(percent_done),'% throguh initial grid search.  Appreciate your patience'))
                    for i6 = 1:numel(mu1_grid)
                        for i7 = 1:numel(mu2_grid)
                            for i8 = 1:numel(mu3_grid)
                                for i9=1:numel(mu4_grid)
                                    for i10=1:numel(mu5_grid)   % dont use n or p as these are carrier concentrations.  q is charge
                                        for i11=1:numel(mu6_grid)   % dont use n or p as these are carrier concentrations.  q is charge
                                            for i12=1:numel(mu7_grid)
                                                for i13=1:numel(mu8_grid)
                                                    current_EF_fixed_mu_vec = [EF_grid(i5) mu1_grid(i6) mu2_grid(i7) mu3_grid(i8) mu4_grid(i9) mu5_grid(i10) mu6_grid(i11) mu7_grid(i12) mu8_grid(i13)];    % create the current EF_fixed_mu_vec
                                                    current_tot_bal = Total_Balance34(current_EF_fixed_mu_vec);
                                                    if current_tot_bal < tot_bal
                                                        EF_fixed_mu_guess34 = current_EF_fixed_mu_vec;
                                                        tot_bal = current_tot_bal;
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end

            case 9
                %         elseif EFFMG34_conditions.num_fixed_elements == 9
                mu1_grid = cell2mat(mu_grid_holder(1));
                mu2_grid = cell2mat(mu_grid_holder(2));
                mu3_grid = cell2mat(mu_grid_holder(3));
                mu4_grid = cell2mat(mu_grid_holder(4));
                mu5_grid = cell2mat(mu_grid_holder(5));
                mu6_grid = cell2mat(mu_grid_holder(6));
                mu7_grid = cell2mat(mu_grid_holder(7));
                mu8_grid = cell2mat(mu_grid_holder(8));
                mu9_grid = cell2mat(mu_grid_holder(9));
                num_grid_points = numel(EF_grid)*numel(mu1_grid)*numel(mu2_grid)*numel(mu3_grid)*numel(mu4_grid)*numel(mu5_grid)*numel(mu6_grid)*numel(mu7_grid)*numel(mu8_grid)*numel(mu9_grid);
                disp(strcat('Starting grid search over_',num2str(num_grid_points),'_points in Ef-mu space'))
                for i5 = 1:numel(EF_grid)
                    percent_done = i5/numel(EF_grid)*100;
                    disp(strcat(num2str(percent_done),'% throguh initial grid search.  Appreciate your patience'))
                    for i6 = 1:numel(mu1_grid)
                        for i7 = 1:numel(mu2_grid)
                            for i8 = 1:numel(mu3_grid)
                                for i9=1:numel(mu4_grid)
                                    for i10=1:numel(mu5_grid)   % dont use n or p as these are carrier concentrations.  q is charge
                                        for i11=1:numel(mu6_grid)   % dont use n or p as these are carrier concentrations.  q is charge
                                            for i12=1:numel(mu7_grid)
                                                for i13=1:numel(mu8_grid)
                                                    for i14=1:numel(mu9_grid)
                                                        current_EF_fixed_mu_vec = [EF_grid(i5) mu1_grid(i6) mu2_grid(i7) mu3_grid(i8) mu4_grid(i9) mu5_grid(i10) mu6_grid(i11) mu7_grid(i12) mu8_grid(i13) mu9_grid(i14) ];    % create the current EF_fixed_mu_vec
                                                        current_tot_bal = Total_Balance34(current_EF_fixed_mu_vec);
                                                        if current_tot_bal < tot_bal
                                                            EF_fixed_mu_guess34 = current_EF_fixed_mu_vec;
                                                            tot_bal = current_tot_bal;
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end

            case 10
                %         elseif EFFMG34_conditions.num_fixed_elements ==
                %         10
                mu1_grid = cell2mat(mu_grid_holder(1));
                mu2_grid = cell2mat(mu_grid_holder(2));
                mu3_grid = cell2mat(mu_grid_holder(3));
                mu4_grid = cell2mat(mu_grid_holder(4));
                mu5_grid = cell2mat(mu_grid_holder(5));
                mu6_grid = cell2mat(mu_grid_holder(6));
                mu7_grid = cell2mat(mu_grid_holder(7));
                mu8_grid = cell2mat(mu_grid_holder(8));
                mu9_grid = cell2mat(mu_grid_holder(9));
                mu10_grid = cell2mat(mu_grid_holder(10));
                num_grid_points = numel(EF_grid)*numel(mu1_grid)*numel(mu2_grid)*numel(mu3_grid)*numel(mu4_grid)*numel(mu5_grid)*numel(mu6_grid)*numel(mu7_grid)*numel(mu8_grid)*numel(mu9_grid)*numel(mu10_grid);
                disp(strcat('Starting grid search over_',num2str(num_grid_points),'_points in Ef-mu space'))
                for i5 = 1:numel(EF_grid)
                    percent_done = i5/numel(EF_grid)*100;
                    disp(strcat(num2str(percent_done),'% throguh initial grid search.  Appreciate your patience'))
                    for i6 = 1:numel(mu1_grid)
                        for i7 = 1:numel(mu2_grid)
                            for i8 = 1:numel(mu3_grid)
                                for i9=1:numel(mu4_grid)
                                    for i10=1:numel(mu5_grid)   % dont use n or p as these are carrier concentrations.  q is charge
                                        for i11=1:numel(mu6_grid)   % dont use n or p as these are carrier concentrations.  q is charge
                                            for i12=1:numel(mu7_grid)
                                                for i13=1:numel(mu8_grid)
                                                    for i14=1:numel(mu9_grid)
                                                        for i15=1:numel(mu10_grid)
                                                            current_EF_fixed_mu_vec = [EF_grid(i5) mu1_grid(i6) mu2_grid(i7) mu3_grid(i8) mu4_grid(i9) mu5_grid(i10) mu6_grid(i11) mu7_grid(i12) mu8_grid(i13) mu9_grid(i14) mu10_grid(i15)];    % create the current EF_fixed_mu_vec
                                                            current_tot_bal = Total_Balance34(current_EF_fixed_mu_vec);

                                                            if current_tot_bal < tot_bal
                                                                EF_fixed_mu_guess34 = current_EF_fixed_mu_vec;
                                                                tot_bal = current_tot_bal;
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end


                % else
            otherwise
                error('looks like more than 10 fixed elements requested - ask developers to add support for more fixed elements in EF_mu_guess34')
                %
                %                 for testing the mu_guess34 procedure -
                %                 plot the errors figure(1) clf
                %                 surf(log10(abs(charge_holder))) figure(2)
                %                 clf surf(log10(abs(element_holder)))
                %                 figure(3) clf
                %                 surf(log10(abs(tot_holder)))

        end  % end the if/elseif  or switch structrue

    end     %%%% end EF_mu_guess for cases 34  %%%%



%% end functions called with dummy names for conditions


end
%% % end main function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
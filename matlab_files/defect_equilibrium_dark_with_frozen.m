% this function calculates the defect equilibrium Energies entered in eV,
% concentrations entered in and come out in #/cm3 Defect charge state
% formation enthalpy is assumed to be in the form dH = dHo + q*EF don't
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
% blocking (1-sum not sum-1 for the finite case).


function [equilib_dark_sol] = defect_equilibrium_dark_with_frozen(conditions, defects)


%%%%%%% build up holder variables for outputs %%%%%%%%%%
N_charge_states_out = zeros(size(conditions.T_equilibrium,2),defects.num_chargestates);   % in the final output, each charge state has a column, each row is a different T. Need this different name for the holder matrix because of scope of n_defects in some of the functions
n_out = zeros(numel(conditions.T_equilibrium),1);   %holder column vector # rows same as T
p_out = zeros(numel(conditions.T_equilibrium),1);
EF_out = zeros(numel(conditions.T_equilibrium),1);
mu_out = zeros(numel(conditions.T_equilibrium),conditions.num_elements);
tot_bal_err_out = zeros(numel(conditions.T_equilibrium),1);
charge_bal_err_out = zeros(numel(conditions.T_equilibrium),1);
element_bal_err_out = zeros(numel(conditions.T_equilibrium),1);
stoich_out = zeros(numel(conditions.T_equilibrium),conditions.num_elements);



%% Tloop_conditions holds *scalar* values for all properties one at a time as the T loop is worked throguh
%%% make a local copy of the conditions variable to use with T dependnet
%%% properties.  At every T we assign a new value of each property to each
%%% variable as needed.  Future fix: assign T dependent dG0 values for
%%% each chargestate for example.
%% all the subroutines will be written to use "conditions_dummy" as a dummy variable but when actually called should be called with "Tloop_conditions"
Tloop_conditions = conditions;

%%%% detect the type of calculation requested.  1=nothing frozen, 2=dfects frozen,
%%%% 3=elements frozen, 4 = defects and elments frozen.  The calc type gets
%%%% written back into the conditions variable
[Tloop_conditions] = Calc_Method(Tloop_conditions);






%%%%%%%%% Main calculation - looped over all the equilibrium T's given.
%%%%%%%%% Tloop_conditions holds only scalar variables and is updated for
%%%%%%%%% each temperature, leaving the original conditions variabl alone
disp('as of now, T dependences for mu and band parameters are being used.  Must add capabilities for other properties being T dependent if desired in future')


%% start the main loop over temperatures
for i1 = 1:numel(conditions.T_equilibrium)  %this instance needs to call back to the original conditions variable that holds the T vector

    Tloop_conditions.T_equilibrium = conditions.T_equilibrium(i1);  % each time through the T loop, set the right values
    Tloop_conditions.kBT_equilibrium = conditions.kBT_equilibrium(i1);
    disp(strcat('Starting calc for T=',num2str(Tloop_conditions.T_equilibrium),'K'))  % call out the start of calc for each T

    % set the mu values for the current temperature.
    Tloop_conditions.mu = conditions.muT_equilibrium(i1,:);  % the mu values set here will be overwritten (ignored) for fixed elements as needed further down in the code

    % Handle set up the ability to specify fixed defect concentrations that
    % vary vs T
    if strcmp(conditions.T_dep_fixed_defect_flag,'On')
        if size(conditions.fixed_defects_concentrations,1)==numel(conditions.T_equilibrium)
            Tloop_conditions.fixed_defects_concentrations = conditions.fixed_defects_concentrations(i1,:);          % take the i1 row of the fixed defects variable from the original conditions file and assign that row in the Tloop_conditions
            Tloop_conditions.fixed_defects_concentrations = Tloop_conditions.fixed_defects_concentrations';  % rotate this to a column vector (doing it in 2 steps to make it obvious we did it)
        else
            error('If T_dep_fixed_defect_flag is ON, then conditions.fixed_defect_concentrations has to be a vector with an entry for each T')
        end
    elseif strcmp(conditions.T_dep_fixed_defect_flag,'Off')
        if size(conditions.fixed_defects_concentrations,2)~=1
            error('If T_dep_fixed_defect_flag is OFF, then conditions.fixed_defect_concentrations has to be a single value (scalar) appropriate for all temperatures')
        end
    end

    % handle T dependent dH values here
    %     elseif strcmp(conditions.T_dep_dH,'On')
    %         Tloop_defects.dH = defects.dHoT(i1,:);
    %
    %     elseif strcmp(conditions.T_dep_dH,'Off')
    %         %% do nothing different
    %     else
    %         error('T dependent fixed defects and or T dependent dHo flags must be on or off')
    %     end
    %

    % Handle T dependent Nc, Nv, and band edges
    if strcmp(conditions.T_dep_bands_flag,'On')
        Tloop_conditions.Eg = conditions.EgT_equilibrium(i1);  % the RHS of these calls need to call into original conditions variable (in case replace all is used by accident)
        Tloop_conditions.Ec = conditions.EcT_equilibrium(i1);
        Tloop_conditions.Ev = conditions.EvT_equilibrium(i1);
        Tloop_conditions.Nc = conditions.NcT_equilibrium(i1);
        Tloop_conditions.Nv = conditions.NvT_equilibrium(i1);
    elseif strcmp(conditions.T_dep_bands_flag,'Off')
        Tloop_conditions.Eg = conditions.EgRef;
        Tloop_conditions.Ec = conditions.EcRef;
        Tloop_conditions.Ev = 0;
        Tloop_conditions.Nc = conditions.NcRef;
        Tloop_conditions.Nv = conditions.NvRef;
    else
        error('T dependent band parameters must be ON or OFF')
    end



    %% run the main calculation according to the method (method won't change vs T)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if Tloop_conditions.calc_method==1 || Tloop_conditions.calc_method ==2  %% These are cases with no elements frozen. 1 is no defects frozen, 2 is some defects frozen.  We get the mu values from the conditions directly


        %% options for the fzero routine - can send these in with conditions variable in future
        % fzero_options = optimset('PlotFcns','optimplotx','MaxFunEvals',1e6,'MaxIter',500,'TolX',1e-4,'TolFun',1e-4);
        % fzero_options = optimset('PlotFcns','optimplotfval','MaxFunEvals',1e6,'MaxIter',5000,'TolX',1e-4,'TolFun',1e-4);
        % fzero_options = optimset('PlotFcns','optimplotfval','MaxFunEvals',1e6,'MaxIter',500,'TolX',5e-3,'TolFun',1e-4);


        disp('Calculation Method 1 or 2 running please be patient....')
        if i1==1
            disp('Starting grid search of Ef to get close to first solution')
            EF_guess = EF_Guess12(Tloop_conditions);   % get a guess for EF close to minimum
            disp('Found guess for first T, passing it off to fzero')
        elseif i1==2
            EF_guess = EF_from_last_T;
            disp('Using solution from prior temperature as guess for solution - make sure to space T values close enough for continuity')
        elseif i1>=3
            % this next line needs to call back to the original conditions
            % variable to have the full list of T values available
            EF_guess = EF_from_last_T + 0.5*(Tloop_conditions.T_equilibrium - last_T)*(EF_from_last_T - EF_from_2nd_last_T)/(last_T - second_last_T);
            disp('Using avg of prior solution and 1st order Taylor expansion as guess...')
        end

        EF_out(i1) = fzero(@Charge_Balance12,EF_guess);
        [n_out(i1), p_out(i1)] = Carrier_Concentrations(EF_out(i1));  %%% carriers
        EF_full_mu_vec = [EF_out(i1) Tloop_conditions.mu];   % create the full EF_mu_vec taking in the fixed mu values
        N_charge_states_out(i1,:) = Chargestate_Concentrations(EF_full_mu_vec);  % compute the chargestate and carrier concentrations from the full EF_mu_vec
        charge_bal_err_out(i1) = Charge_Balance12(EF_out(i1));   % compute the net charge just to check final answers
        mu_out(i1,:) = Tloop_conditions.mu;

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


    elseif Tloop_conditions.calc_method==3 || Tloop_conditions.calc_method==4 %% These are the cases where some elements are frozen.  3 is only elements frozen, 4 is elements and defects frozen
        disp('Calculation Method 3 or 4 running please be patient....')
        disp('WARNING!!!: Remember, when element concentrations are fixed, it is possible that the desired solution is impossible to find within the limits of mu the user specified for each fixed element.  So check the results carefully and spend some time bracketing the solution by hand.  Remember that only complexes couple impurity elements together so you can vary one mu at a time for the most part. ')

        if strcmp(conditions.search_method_flag,'particleswarm_pattern') || strcmp(conditions.search_method_flag,'particleswarm_pattern_simplex')
            nvars = Tloop_conditions.num_fixed_elements + 1;
            iter_lim = Tloop_conditions.fixed_element_swarm_iterlimit_base*nvars;
            stall_lim = Tloop_conditions.fixed_element_swarm_stall_lim;
            
            if i1 == 1
                EF_min = Tloop_conditions.Ev - 5*Tloop_conditions.kBT_equilibrium;
                EF_max = Tloop_conditions.Ec + 5*Tloop_conditions.kBT_equilibrium;
                lb = [EF_min; Tloop_conditions.fixed_elements_mu_ranges(Tloop_conditions.indices_of_fixed_elements,1)];
                ub = [EF_max; Tloop_conditions.fixed_elements_mu_ranges(Tloop_conditions.indices_of_fixed_elements,2)];
                particles_per_fixed_mu = (ub-lb) / Tloop_conditions.kBT_equilibrium;
                grid_product = prod(particles_per_fixed_mu);  % number of points needed to fill the search hypervolume spaced kBT apart
                swarmsize1 = ceil(min(Tloop_conditions.fixed_element_swarm1_max_size, grid_product^((nvars-1)/nvars))); % initiate a very large swarm to densely sample the space.  the idea is that you want to scale better than the dimensionality of the space being searched
                Tloop_conditions.ps_options1 = optimoptions('particleswarm','SwarmSize',swarmsize1,'Display','iter','ObjectiveLimit',1e3,'MaxIterations',iter_lim,'MaxStallIterations',stall_lim,'HybridFcn', @patternsearch);
                tic
                EF_fixed_mu_vec_sol = particleswarm(@Total_Balance34,nvars,lb,ub,Tloop_conditions.ps_options1);
                if strcmp(conditions.search_method_flag,'particleswarm_pattern_simplex')
                    fminsearch_options = optimset('MaxFunEvals',10000,'MaxIter',5000,'TolX',1e-6,'TolFun',1e-6);  % may need to tune these settings
                    EF_fixed_mu_vec_sol = fminsearch(@Total_Balance34,EF_fixed_mu_vec_sol, fminsearch_options);
                end
                guess_time = toc;
            elseif i1>=2
                lb = EF_fixed_mu_from_last_T - Tloop_conditions.fixed_element_swarm2_search_band_kB*Tloop_conditions.kBT_equilibrium;
                ub = EF_fixed_mu_from_last_T + Tloop_conditions.fixed_element_swarm2_search_band_kB*Tloop_conditions.kBT_equilibrium;
                particles_per_fixed_mu = 2*Tloop_conditions.fixed_element_swarm2_search_band_kB;
                grid_product = prod(particles_per_fixed_mu);  % number of points needed to fill the search hypervolume spaced kBT apart
                swarmsize2 = ceil(max(Tloop_conditions.fixed_element_swarm2_min_size, grid_product^((nvars-1)/nvars))); % initiate a very large swarm to densely sample the space.  the idea is that you want to scale better than the dimensionality of the space being searched
                Tloop_conditions.ps_options2 = optimoptions('particleswarm','SwarmSize',swarmsize2,'Display','iter','ObjectiveLimit',1e3,'MaxIterations',iter_lim,'MaxStallIterations',stall_lim,'HybridFcn', @patternsearch);
                tic
                EF_fixed_mu_vec_sol = particleswarm(@Total_Balance34,nvars,lb,ub,Tloop_conditions.ps_options2);
                if strcmp(conditions.search_method_flag,'particleswarm_pattern_simplex')
                    fminsearch_options = optimset('MaxFunEvals',10000,'MaxIter',500,'TolX',1e-6,'TolFun',1e-6);  % may need to tune these settings
                    EF_fixed_mu_vec_sol = fminsearch(@Total_Balance34,EF_fixed_mu_vec_sol, fminsearch_options);
                end
                guess_time = toc;
            end

        elseif strcmp(conditions.search_method_flag,'grid_fminsearch')
            max_mu_range = max(abs(Tloop_conditions.fixed_elements_mu_ranges(:,2)-Tloop_conditions.fixed_elements_mu_ranges(:,1)));
            start_kBT = conditions.T_equilibrium(1)*Tloop_conditions.kB;  %waht is kBT at the first temperature where we do the grid search?
            fine_factor = 8;
            Tloop_conditions.fixed_elements_mu_npoints = ceil(fine_factor*max_mu_range/start_kBT);
            fminsearch_options = optimset('MaxFunEvals',Tloop_conditions.fixed_element_fmin_MaxFunEvals,'MaxIter',Tloop_conditions.fixed_element_fmin_MaxIter,'TolX',Tloop_conditions.fixed_element_fmin_TolX,'TolFun',Tloop_conditions.fixed_element_fmin_TolFun);  % may need to tune these settings
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
                EF_fixed_mu_guess = EF_fixed_mu_from_last_T + 0.5*((EF_fixed_mu_from_last_T - EF_fixed_mu_from_2nd_last_T)/(conditions.T_equilibrium(i1-1) - conditions.T_equilibrium(i1-2)))*(conditions.T_equilibrium(i1) - conditions.T_equilibrium(i1-1));
                disp('Using avg of prior solution and 1st order Taylor expansion as guess...')
            end

            % do the optimization using fminsearch
            tic
            [EF_fixed_mu_vec_sol,~,fminsearch_exit_flag] = fminsearch(@Total_Balance34,EF_fixed_mu_guess, fminsearch_options);  % up to here we deal with EF_fixed_mu_vec.  After the optimization we rebuild the whole EF_full_mu_vec
            [EF_fixed_mu_vec_sol,~,fminsearch_exit_flag] = fminsearch(@Total_Balance34,EF_fixed_mu_vec_sol, fminsearch_options);
            toc


        else
            error('For case 3 & 4 calculations, search_method_flag must be one of the coded options')
        end


        % these items here are executed for cases 3&4 using either solver option
        disp(strcat('Found solution as__',num2str(EF_fixed_mu_vec_sol),'__for EF_fixed_mu_vec'))
        EF_full_mu_vec_sol = Expand_Fixed_Mu_Vec_To_Full_Mu_Vec(EF_fixed_mu_vec_sol);  % expand it out to the full EF_mu_vec
        EF_out(i1) = EF_full_mu_vec_sol(1);  % pick off the optimized EF value
        mu_out(i1,:) = EF_full_mu_vec_sol(2:end);
        [n_out(i1), p_out(i1)]  = Carrier_Concentrations(EF_out(i1));  %%% carriers
        N_charge_states_out(i1,:) = Chargestate_Concentrations(EF_full_mu_vec_sol);  %% chargestates
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

% there was a little loop here indexed by i2 but deleted it
%%%%%%%%%%%%%% end of main calc %%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%% create final output strcutrues  %%%%%%%%%%%%%%%%%%%%%%%%

%% sum over chargstates in each defect to get the totals for each defect. Do this here at the end so the loop only runs one time %%
for i3=1:defects.num_defects
    index = defects.cs_ID==i3;
    N_defects_out(:,i3) = sum(N_charge_states_out(:,index),2);
end


%% build up the solution structure
equilib_dark_sol.defect_names = defects.defect_names;
equilib_dark_sol.chargestate_names = defects.chargestate_names;
equilib_dark_sol.T_equilibrium = (conditions.T_equilibrium)';
equilib_dark_sol.Nd = conditions.Nd*ones(numel(conditions.T_equilibrium),1);
equilib_dark_sol.Na = conditions.Na*ones(numel(conditions.T_equilibrium),1);
equilib_dark_sol.EFn = EF_out;   % planning ahead for light calcs
equilib_dark_sol.EFp = EF_out;
equilib_dark_sol.n = n_out;
equilib_dark_sol.p = p_out;
equilib_dark_sol.chargestates = N_charge_states_out;
equilib_dark_sol.defects = N_defects_out;
[equilib_dark_sol.Ga_O_stoich, equilib_dark_sol.elements(:,1), equilib_dark_sol.elements(:,2), equilib_dark_sol.elements(:,3), equilib_dark_sol.elements(:,4), equilib_dark_sol.elements(:,5), equilib_dark_sol.elements(:,6), equilib_dark_sol.elements(:,7), equilib_dark_sol.elements(:,8), equilib_dark_sol.elements(:,9), equilib_dark_sol.elements(:,10), equilib_dark_sol.elements(:,11), equilib_dark_sol.elements(:,12), equilib_dark_sol.elements(:,13), equilib_dark_sol.elements(:,14), equilib_dark_sol.elements(:,15), equilib_dark_sol.elements(:,16), equilib_dark_sol.elements(:,17)] =  Ga2O3_stoich(equilib_dark_sol, conditions, defects);
equilib_dark_sol.charge_bal_err = charge_bal_err_out;
equilib_dark_sol.element_bal_err = element_bal_err_out;
equilib_dark_sol.tot_bal_err = tot_bal_err_out;
equilib_dark_sol.mu = mu_out;


%%%%%%%%% end of main function %%%%%%%%%%%%%%%




%% nested subroutines - note these all share memory space with the main function ("end" for the main function comes after all of these subroutines are ended %%

%% These functions have the conditions variable as input.  When actually calling these functions we want to call it with "Tloop_conditions" as input.  They are written with a dummy variable "local conditions" to differentiate the context

    function [conditions_dummy] = Calc_Method(conditions_dummy)
        %% check for valid imputs for fixed and open defects and elements
        if max(size(conditions_dummy.fixed_defects))~=defects.num_defects  % check if each defect has a 1 or 0
            error('Tloop_conditions.fixed_defects has to be same size as num_defects')
        elseif max(size(conditions_dummy.fixed_elements))~=defects.numelements  % check that each element part of the calculation is listed
            error('Tloop_conditions.fixed_elements has to be same size as # elements')
        elseif sum(conditions_dummy.fixed_defects) ~= conditions_dummy.num_fixed_defects
            error('Something wrong with number of fixed defects')
        elseif sum(conditions_dummy.fixed_elements) ~= conditions_dummy.num_fixed_elements
            error('Something wrong with number of fixed elements')
        end

        %% these logic checks figure out how calc will be done
        if sum(conditions_dummy.fixed_defects)==0   % no defects fixed, all open
            conditions_dummy.some_defects_fixed_flag = 0;
            conditions_dummy.fixed_defects_index = [];
            conditions_dummy.open_defects_index = ones(size(conditions_dummy.fixed_defects));
            conditions_dummy.num_fixed_defects = 0;
            conditions_dummy.num_open_defects = max(size(conditions_dummy.open_defects_index));
        elseif sum(conditions_dummy.fixed_defects)~=0  %% some defects have fixed concentration
            conditions_dummy.some_defects_fixed_flag = 1;
            conditions_dummy.fixed_defects_index = find(conditions_dummy.fixed_defects==1);
            conditions_dummy.open_defects_index = find(conditions_dummy.fixed_defects==0);
            conditions_dummy.num_fixed_defects = max(size(conditions_dummy.fixed_defects_index));
            conditions_dummy.num_open_defects = max(size(conditions_dummy.open_defects_index));
        end

        if sum(conditions_dummy.fixed_elements)==0   % no elements fixed, all open (mu set for all not concentration)
            conditions_dummy.some_elements_fixed_flag = 0;
            conditions_dummy.fixed_elements_index = [];
            conditions_dummy.open_elements_index = ones(size(conditions_dummy.fixed_elements));
        elseif sum(conditions_dummy.fixed_elements)~=0  %% at least one element concentration is set
            conditions_dummy.some_elements_fixed_flag = 1;
            conditions_dummy.fixed_elements_index = find(conditions_dummy.fixed_elements==1);
            conditions_dummy.open_elements_index = find(conditions_dummy.fixed_elements==0);
            conditions_dummy.num_fixed_elements = numel(conditions_dummy.fixed_elements_index);
        end

        if conditions_dummy.some_defects_fixed_flag==0 && conditions_dummy.some_elements_fixed_flag==0
            conditions_dummy.calc_method = 1;
        elseif conditions_dummy.some_defects_fixed_flag==1 && conditions_dummy.some_elements_fixed_flag==0
            conditions_dummy.calc_method = 2;
        elseif conditions_dummy.some_defects_fixed_flag==0 && conditions_dummy.some_elements_fixed_flag==1
            conditions_dummy.calc_method = 3;
        elseif conditions_dummy.some_defects_fixed_flag==1 && conditions_dummy.some_elements_fixed_flag==1
            conditions_dummy.calc_method = 4;
        else
            error('defects and elements must be either fixed or open')
        end

    end



    function [EF_guess12] = EF_Guess12(conditions_dummy)
        %%%% this will only be called for mthods 1 & 2. for methods 1&2,
        %%%% function that uses grid search to get an EF guess close to the
        %%%% charge balance solution.  It should return 1 value (old method
        %%%% tried to give 2 but it was overcomplicated - leave it as
        %%%% comments for now in case it is needed sometime)
        % keep kBT and Eg as local variables in this function by using
        % distinct aliases for the input arguments

        [EF_grid] = EF_Grid_Maker(conditions_dummy);

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



    function [EF_fixed_mu_guess34] = EF_Fixed_Mu_Guess34(conditions_dummy)
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
        [EF_grid] = EF_Grid_Maker(conditions_dummy);

        % create search grid for each frozen element
        [mu_grid_holder] = Fixed_Mu_Ggrid_Maker(conditions_dummy);   % this creates a cell array with number of cells equal to the number of fixed elements.  In those cells are vectors of mu grids for each fixed element.

        % start going through the cases for diffrent numbers of fixed
        % elements
        tot_bal = Inf;  %initialize tot_bal as infinite (< operator will pick any finite value as less than on first iteration)

        % converting this to a switch structrue so it only evaluates one
        % case rather than going through all the elseifs, since this gets
        % called a lot


        switch conditions_dummy.num_fixed_elements
            case 1

                %         if conditions_dummy.num_fixed_elements == 1
                mu1_grid = cell2mat(mu_grid_holder(1));   % extract the first mu vec to search over and make it a vector
                %                 diagnstics
                %                 charge_holder = zeros(numel(EF_grid),numel(mu1_grid));
                %                 element_holder = zeros(numel(EF_grid),numel(mu1_grid));
                %                 tot_holder = zeros(numel(EF_grid),numel(mu1_grid));

                num_grid_points = numel(EF_grid)*numel(mu1_grid);
                disp(strcat('Starting grid search over_',num2str(num_grid_points),'_points in Ef-mu space'))
                for i5 = 1:numel(EF_grid)
                    for i6 = 1:numel(mu1_grid)
                        current_EF_fixed_mu_vec = [EF_grid(i5) mu1_grid(i6)];         % create the current EF_fixed_mu_vec
                        current_tot_bal = Total_Balance34(current_EF_fixed_mu_vec); % throw it into the total balance function
                        %                         charge_holder(i5,i6) = Charge_Balance34(current_EF_fixed_mu_vec);
                        %                         element_holder(i5,i6) = sum(Fixed_Element_Balance34(current_EF_fixed_mu_vec));
                        %                         tot_holder(i5,i6) = abs(charge_holder(i5,i6)) + abs(element_holder(i5,i6));
                        if current_tot_bal < tot_bal    % see if it is better than prior values.  If yes, set guess to this
                            EF_fixed_mu_guess34 = current_EF_fixed_mu_vec;
                            tot_bal = current_tot_bal;
                        end
                    end
                end

            case 2
                %         elseif conditions_dummy.num_fixed_elements == 2
                mu1_grid = cell2mat(mu_grid_holder(1));
                mu2_grid = cell2mat(mu_grid_holder(2));
                num_grid_points = numel(EF_grid)*numel(mu1_grid)*numel(mu2_grid);
                %                 update_progress_num = floor(num_grid_points/10);
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
                %         elseif conditions_dummy.num_fixed_elements == 3
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
                %         elseif conditions_dummy.num_fixed_elements == 4
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
                %         elseif conditions_dummy.num_fixed_elements == 5
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
                %         elseif conditions_dummy.num_fixed_elements == 6
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
                %         elseif conditions_dummy.num_fixed_elements == 7
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
                %         elseif conditions_dummy.num_fixed_elements == 8
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
                %         elseif conditions_dummy.num_fixed_elements == 9
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
                %         elseif conditions_dummy.num_fixed_elements == 10
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
                %                 for testing the mu_guess34 procedure - plot the errors
                %                 figure(1)
                %                 clf
                %                 surf(log10(abs(charge_holder)))
                %                 figure(2)
                %                 clf
                %                 surf(log10(abs(element_holder)))
                %                 figure(3)
                %                 clf
                %                 surf(log10(abs(tot_holder)))

        end  % end the if/elseif  or switch structrue

    end     %%%% end EF_mu_guess for cases 34  %%%%



    function [fixed_mu_grids]  = Fixed_Mu_Ggrid_Maker(conditions_dummy)
        % this gets called only for methods 3 &4, and gets called before
        % anything else function to create an array holding vectors of mu
        % values for each frozen element (and no others).  The fixed
        % elements will be in order in the array, but not the same absolute
        % element numbers. So if there are 6 elements, and 1, 3, and 5 are
        % frozen, this will have 3 cells, 1=1st, 3=2nd , 5=3rd
        %% future possible improvement: can just make holder a matrix directly rather than making a cell array here and turning it back into a vector when used in guess34
        %% also could move this task into the script instead and include this object in the conditions variable sent.  In that case it would be better to have it as a cell array as then some entries can be nulls rathet than zeros required if this is a matrix

        fixed_mu_grids = cell(1,conditions_dummy.num_fixed_elements);  % create empty cell array with enough entries to hold a vector of mu values for each fixed element
        counter = 1;
        for i16 = conditions_dummy.indices_of_fixed_elements  % loop over the i16 values indexing the frozen elements  This way only executes the loop for the frozen ones, skipping the non-frozen ones so gets rid of need to insert if elseif inside the for loop.  This will only be called in methods 3 and 4 so no need to worry about methods 1&2 when none are frozen
            mu_inc = (conditions_dummy.fixed_elements_mu_ranges(i16,2)-conditions_dummy.fixed_elements_mu_ranges(i16,1))/conditions_dummy.fixed_elements_mu_npoints; % decide up front that you will divide the rnage into a fixed number of increments
            fixed_mu_grids{counter} = conditions_dummy.fixed_elements_mu_ranges(i16,1) : mu_inc : conditions_dummy.fixed_elements_mu_ranges(i16,2);
            counter = counter + 1;
        end
    end



    function [EF_grid] = EF_Grid_Maker(conditions_dummy)
        % Set up a grid of EF values the grid.  Using interval 1/2 of kBT
        % guarantees you can't miss the solution whcih should be thus
        % within kB/2 of the guess
        EF_inc = conditions_dummy.kBT_equilibrium/2;
        EF_grid = conditions_dummy.Ev-10*EF_inc : EF_inc : ceil((conditions_dummy.Ec-conditions_dummy.Ev)/EF_inc)*EF_inc + 10*EF_inc;  % this makes a grid to check over the range -5kBT to Eg+5kBT, the ceil() command takes care of casee where Ec-Ev/interval is not an integer
    end

%% end functions called with Tloop_conditions



%% These functions are involved in calculating the numbers of things.  Some take in reduced and others full EF_mu_vecs holding just the fixed mu's or the full list of mu's respectively

    function [charge_bal12] = Charge_Balance12(EF_dummy)
        % function to calculate signed charge balance for methods 1 and 2
        % for speedup.  Same as the 34 version but stripped all the fixed
        % mu stuff out.  This gets called a lot so dont put any logic
        % checks inside - do checking outside before calling this
        [n,p] = Carrier_Concentrations(EF_dummy);  %deal with carriers first
        EF_full_mu_vec_dummy = [EF_dummy Tloop_conditions.mu];   % create the full EF_mu_vec
        [N_chargestates] = Chargestate_Concentrations(EF_full_mu_vec_dummy);
        charge_bal12 = sum(defects.cs_charge.* N_chargestates') + p - n + Tloop_conditions.Nd - Tloop_conditions.Na  ;              % absolute signed charge bal comes out signed +/-
        % charge_bal12 = (sum(defects.cs_charge.* N_chargestates') + p - n + Tloop_conditions.Nd - Tloop_conditions.Na)/sqrt(Tloop_conditions.Nc*Tloop_conditions.Nv)  ;  % scaled to ni2 sort of.   comes out signed +/-
    end %end charge_bal12



    function [charge_bal34] = Charge_Balance34(EF_fixed_mu_vec_dummy)
        % function to calculate signed charge balance.  This gets called a
        % lot so dont put any logic checks inside - do checking outside
        % before calling this
        [n,p] = Carrier_Concentrations(EF_fixed_mu_vec_dummy(1));  %deal with carriers first
        % now deal with defects.  Need to take in reduced mu vec and expand
        % out to the full one before sending to chargestate calculation
        [EF_full_mu_vec_dummy] = Expand_Fixed_Mu_Vec_To_Full_Mu_Vec(EF_fixed_mu_vec_dummy);
        [N_chargestates] = Chargestate_Concentrations(EF_full_mu_vec_dummy);  % send this full vec to calculate the chargestates
        charge_bal34 = sum(defects.cs_charge.* N_chargestates') + p - n + Tloop_conditions.Nd - Tloop_conditions.Na;   % absolute signed charge bal
        % charge_bal34 = (sum(defects.cs_charge.* N_chargestates') + p - n + Tloop_conditions.Nd - Tloop_conditions.Na)/sqrt(Tloop_conditions.Nc*Tloop_conditions.Nv);  % relative comes out signed +/-
    end %end charge_bal34



    function [n,p] = Carrier_Concentrations(EF_dummy)
        %%%% function to compute n and p.  Input must be EF only (scalar
        %%%% not vector with mu's attached)  %%%%%%%%%%%%
        etaCB = (EF_dummy - Tloop_conditions.Ec)/Tloop_conditions.kBT_equilibrium;
        etaVB = (EF_dummy - Tloop_conditions.Ev)/Tloop_conditions.kBT_equilibrium;   % this looks wrong (not symmetric compared to CB case) but it is right. Direction of integration and sign on EF-Ev are both swapped.

        if strcmp(Tloop_conditions.Boltz_or_FD_flag,'FD')            % use Fermi-Dirac integrals so degenerate conditions handled correctly
            [n,~] = n_Fermi_Dirac(etaCB,Tloop_conditions.Nc);
            [p,~] = n_Fermi_Dirac(-etaVB,Tloop_conditions.Nv);
        elseif strcmp(Tloop_conditions.Boltz_or_FD_flag,'Boltz')              % %     use just Boltzmann approx
            n = Tloop_conditions.Nc*exp(etaCB);   % these are right (sign swap).  Boltzmann factors should end up <1 when EF is in gap
            p = Tloop_conditions.Nv*exp(-etaVB);
        else
            error('Boltz_or_FD_flag must be Boltz or FD')
        end
    end



    function [tot_bal34] = Total_Balance34(EF_fixed_mu_vec_dummy)
        %% this constructs the objective function that is to be zeroed out.  Takes in the reduced EF_mu_Vec and gives a scalar out
        % This is the log10 of the total error from zero on the element sum
        % and charge sum
        %         tot_bal34 = log10(abs(Charge_Balance34(EF_fixed_mu_vec_dummy)) + sum(abs(Fixed_Element_Balance34(EF_fixed_mu_vec_dummy)),'omitnan'));
        %         tot_bal34 = log10(abs(Charge_Balance34(EF_fixed_mu_vec_dummy)) + sum(abs(Fixed_Element_Balance34(EF_fixed_mu_vec_dummy))));
        tot_bal34 = abs(Charge_Balance34(EF_fixed_mu_vec_dummy)) + sum(abs(Fixed_Element_Balance34(EF_fixed_mu_vec_dummy)));
    end  %%% end tot_balance objective function  %%%%%%%%%%%%%




    function [fixed_element_bal34_vec] = Fixed_Element_Balance34(EF_fixed_mu_vec_dummy)
        % returns a vector of the signed difference between fixed element
        % concs and their targets, only for the fixed elements
        %         element_bal34 = zeros(size(Tloop_conditions.num_fixed_elements));  % initialize as all zeros
        [EF_full_mu_vec_dummy] = Expand_Fixed_Mu_Vec_To_Full_Mu_Vec(EF_fixed_mu_vec_dummy);
        [N_chargestates] = Chargestate_Concentrations(EF_full_mu_vec_dummy);  % send this full vec to calculate the chargestates
        %         fixed_element_bal34_vec =  ((((defects.cs_dm')*N_chargestates')' - Tloop_conditions.fixed_elements_concentrations).*Tloop_conditions.fixed_elements)./Tloop_conditions.fixed_elements_concentrations;   % have to rotate and reshape the arrays a few times to get sizes and shapes right.  Then .* by the fixed element mask to isolate only the target elements - for example a substitutional defect also produces a missing host atom so dont want to include those. This masking is applied
        fixed_element_bal34_vec =  (((defects.cs_dm')*N_chargestates')' - Tloop_conditions.fixed_elements_concentrations).*Tloop_conditions.fixed_elements;   % have to rotate and reshape the arrays a few times to get sizes and shapes right.  Then .* by the fixed element mask to isolate only the target elements - for example a substitutional defect also produces a missing host atom so dont want to include those. This masking is applied
    end % end element bal34



    function [N_chargestates] = Chargestate_Concentrations(EF_full_mu_vec_dummy)
        %%%%  function that calculates the numbers of each charge state of
        %%%%  each defect.  Input has to be a full EF_mu_vector with
        %%%%  entries for ALL elements not just fixed ones.  Fixed defects
        %%%%  overrides fixed elements, thus the 4 method cases are spit
        %%%%  according to whether or not any defects are fixed.  The main
        %%%%  calc loop takes care of supplying a full EF_mu_vec to send
        %%%%  in.

        if Tloop_conditions.calc_method==1 || Tloop_conditions.calc_method==3

            if strcmp(Tloop_conditions.vib_ent_flag,'On')
                dG = defects.cs_dHo - defects.cs_dm*EF_full_mu_vec_dummy(2:end)' + defects.cs_charge*EF_full_mu_vec_dummy(1) - 3*Tloop_conditions.kBT_equilibrium*sum(defects.cs_dm,2);  % last term is -TdS.  A vacacny has sum(cs_dm)=-1, and interstitial has +1 so dS is + for interstitials
            elseif strcmp(Tloop_conditions.vib_ent_flag,'Off')
                dG = defects.cs_dHo - defects.cs_dm*EF_full_mu_vec_dummy(2:end)' + defects.cs_charge*EF_full_mu_vec_dummy(1);
            else
                error('Vibrational entropy flag must be On or Off')
            end

            if strcmp(Tloop_conditions.site_blocking_flag,'On_Infinite')
                N_chargestates = (defects.cs_prefactor .* exp(-dG/Tloop_conditions.kBT_equilibrium) ./ (1+sum(exp(-dG/Tloop_conditions.kBT_equilibrium)) ))';    %% with site blocking for infinite lattice
            elseif strcmp(Tloop_conditions.site_blocking_flag,'On_Finite')
                N_chargestates = (defects.cs_prefactor .* exp(-dG/Tloop_conditions.kBT_equilibrium) ./ (1-sum(exp(-dG/Tloop_conditions.kBT_equilibrium)) ))';    %% with site blocking for finite crystal
            elseif strcmp(Tloop_conditions.site_blocking_flag,'Off')
                N_chargestates = (defects.cs_prefactor .* exp(-dG/Tloop_conditions.kBT_equilibrium))';   % Boltzmann dilute defects approximation here for all charge states of all defects.
            else
                error('Site blocking flag must be On_Infinite, On_Finite, or Off')
            end


        elseif Tloop_conditions.calc_method==2 || Tloop_conditions.calc_method==4
            % total formation energy with chemical potential part
            if strcmp(Tloop_conditions.vib_ent_flag,'On')
                dG = defects.cs_dHo - defects.cs_dm*EF_full_mu_vec_dummy(2:end)' + defects.cs_charge*EF_full_mu_vec_dummy(1) - 3*Tloop_conditions.kBT_equilibrium*sum(defects.cs_dm,2);  % last term is -TdS.  A vacacny has sum(cs_dm)=-1, and interstitial has +1 so dS is + for interstitials
            elseif strcmp(Tloop_conditions.vib_ent_flag,'Off')
                dG = defects.cs_dHo - defects.cs_dm*EF_full_mu_vec_dummy(2:end)' + defects.cs_charge*EF_full_mu_vec_dummy(1);
            else
                error('Vibrational entropy flag must be On or Off')
            end

            % relative formation energy without chemical potential part. We
            % will only use it for the fixed ones
            N_chargestates = zeros(size(dG));  % initialize an empty holder because we will take care of the fixed ones first
            if strcmp(Tloop_conditions.vib_ent_flag,'On')
                dG_rel = defects.cs_dHo + defects.cs_charge*EF_full_mu_vec_dummy(1) - 3*Tloop_conditions.kBT_equilibrium*sum(defects.cs_dm,2);  % last term is -TdS.  A vacacny has sum(cs_dm)=-1, and interstitial has +1 so dS is + for interstitials
            elseif strcmp(Tloop_conditions.vib_ent_flag,'Off')
                dG_rel = defects.cs_dHo + defects.cs_charge*EF_full_mu_vec_dummy(1);
            else
                error('Vibrational entropy flag must be On or Off')
            end

            Boltz_fac_rel = exp(-dG_rel/Tloop_conditions.kBT_equilibrium); % Boltzmann factors for all chargstates - to be used in Gibbs distributions for each fixed/frozen defects

            % Here we create effective dG values for the fixed defects and
            % overwrite the values for the open case but just for the fixed
            % defects
            for i17 = 1:Tloop_conditions.num_fixed_defects     % loop through the fixed defects one at a time and assign their chargestates first
                which_frozen_defects_index = Tloop_conditions.fixed_defects_index(i17);
                which_frozen_chargestates_index = defects.cs_indices_lo(which_frozen_defects_index):defects.cs_indices_hi(which_frozen_defects_index);
                Z = sum(Boltz_fac_rel(which_frozen_chargestates_index));  % this is the Z for this particular i17'th defect
                N_chargestates(which_frozen_chargestates_index) = Tloop_conditions.fixed_defects_concentrations(which_frozen_defects_index)*Boltz_fac_rel(which_frozen_chargestates_index)/Z ;     % denom is a scalar.  Compute the conc of each charge state in defect i6 is Boltz factor/Z(i)*total conc  Totalconc is a scalar also - same for all chargestates of the give defetc
                dG(which_frozen_chargestates_index) = -Tloop_conditions.kBT_equilibrium*log(N_chargestates(which_frozen_chargestates_index)./defects.cs_prefactor(which_frozen_chargestates_index));   % compute the dG that would give the concentration we just calculated and repalce the one calculated the normal way
            end

            % at this point , we have dG values for all chargestates that
            % will give the right concentrations - pass them on now to the
            % normal calc routine to do all of them at once.

            if strcmp(Tloop_conditions.site_blocking_flag,'On_Infinite')
                N_chargestates = (defects.cs_prefactor .* exp(-dG/Tloop_conditions.kBT_equilibrium) ./ (1+sum(exp(-dG/Tloop_conditions.kBT_equilibrium)) ))';    %% with site blocking for infinite lattice
            elseif strcmp(Tloop_conditions.site_blocking_flag,'On_Finite')
                N_chargestates = (defects.cs_prefactor .* exp(-dG/Tloop_conditions.kBT_equilibrium) ./ (1-sum(exp(-dG/Tloop_conditions.kBT_equilibrium)) ))';    %% with site blocking for finite crystal
            elseif strcmp(Tloop_conditions.site_blocking_flag,'Off')
                N_chargestates = (defects.cs_prefactor .* exp(-dG/Tloop_conditions.kBT_equilibrium))';   % Boltzmann dilute defects approximation here for all charge states of all defects.
            else
                error('Site blocking flag must be On_Infinite, On_Finite, or Off')
            end
        else
            error('must be one of four calculation methods based on fixed/open elements and defects')
        end
    end   % end chargestate_concentrations



    function [EF_full_mu_vec_out] = Expand_Fixed_Mu_Vec_To_Full_Mu_Vec(EF_fixed_mu_vec_dummy)
        % does what the name says - takes in a reduced EF_fixed_mu_vec and
        % spits out the whole one with mu values for all elements.
        fixed_mu_vec = EF_fixed_mu_vec_dummy(2:end);  %strip off EF from front so only have the fixed mu's
        full_mu_vec = Tloop_conditions.mu;   % initialize as those currently in the Tloop_conditions.mu - this assigns all the NON-fixed ones as well as the fixed ones
        counter = 0;
        for i18 = Tloop_conditions.indices_of_fixed_elements  % loop calling only the fixed elements
            counter = counter + 1;
            full_mu_vec(i18) = fixed_mu_vec(counter);  % write fixed mu's into the full mu vec.
        end
        EF_full_mu_vec_out = [EF_fixed_mu_vec_dummy(1) full_mu_vec];
    end


end  %%% end whole main function


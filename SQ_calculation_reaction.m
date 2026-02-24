function [SQ_dark_sol, SQ_conditions] = SQ_calculation_reaction(dummy_defects, dummy_conditions)
% important to do: need to make it so it calcs the mu values at the
% current_Tfreeze for T dependent mu's or those set by 2nd phases

figure(1)
hold on

% check to make sure the hi, lo, and increments are decresing
if dummy_conditions.SQ_T_step>0
    dummy_conditions.SQ_T_step = -dummy_conditions.SQ_T_step;
elseif dummy_conditions.SQ_T_start < dummy_conditions.SQ_T_end  % flip the start and end values
    temp = dummy_conditions.SQ_T_end;
    dummy_conditions.SQ_T_end = dummy_conditions.SQ_T_start;
    dummy_conditions.SQ_T_start = temp;
    clear temp
elseif dummy_conditions.SQ_T_end==dummy_conditions.SQ_T_start
    error('Tstart and Tend can not be the same')
elseif dummy_conditions.SQ_T_step == 0
    error("SQ_T_step can't be 0")
end

%% working and output conditions variables
% take in dummy_conditions, we need to keep it throughout the whole calc
% becasue we end up sending working_conditions to the calc engine with single T
% values at a time.  So, as we go down in T for each i,j we need to reset all the defects to their original frozen/open state.
% So dont delete or modify the incoming dummy_conditions.
% also copy dummy_conditions to SQ_conditions that we will output
% so yes it makes sense to make 2 copies of the variable here.
SQ_conditions = dummy_conditions;
working_conditions = dummy_conditions;

% erase the Tequilibrium and T dependent properties as we will replace them
% with the list of values at the actual temperatures we calculate at.  The
% dummy_conditions contains the 1 K spaced values we will pull from
SQ_conditions.T_equilibrium = [];
SQ_conditions.num_T_equilibrium =[];
SQ_conditions.kBT_equilibrium = [];
SQ_conditions.muT_equilibrium = [];
SQ_conditions.NcT_equilibrium = [];
SQ_conditions.NvT_equilibrium = [];
SQ_conditions.EgT_equilibrium = [];
SQ_conditions.EcT_equilibrium = [];
SQ_conditions.EvT_equilibrium = [];

working_conditions.T_equilibrium = [];
working_conditions.num_T_equilibrium =[];
working_conditions.kBT_equilibrium = [];
working_conditions.muT_equilibrium = [];
working_conditions.NcT_equilibrium = [];
working_conditions.NvT_equilibrium = [];
working_conditions.EgT_equilibrium = [];
working_conditions.EcT_equilibrium = [];
working_conditions.EvT_equilibrium = [];
%% end initializing two copies of conditions variables


% in the script that calls this pogram, we define 1 K steps from Thigh to Tlow so that we get all
% the needed things like mu and Eg for every T in 1 K spaced intervals.
% But we dont necessarily want to do the full SQ calc starting from each 1
% K - probably we want it at like every 100 K or 50 K.  So these
% temperatures give the ones we will atually start from.
% these two things dont exist in the incoming conditions vector - create new
% variables for the list of temperatures we will calculate at evenly spaced
% from Thi to Tlow

T_reg_spaced = dummy_conditions.SQ_T_start : dummy_conditions.SQ_T_step : dummy_conditions.SQ_T_end;
num_T_reg_space = numel(T_reg_spaced);

% the max possible number of T values we will end up calculating at could
% be this big.  Each i,j pair will probably end up having fewer though.
max_possible_SQ_T_values = dummy_defects.num_defects + num_T_reg_space;

% for each ij, each defect will end up with one Tfreeze value
Tfreeze_ij_for_defect_k_holder = zeros(dummy_conditions.SQ_num_chardists, dummy_conditions.SQ_num_Trates, dummy_defects.num_defects);

% initialize holder variables for things we will write to both SQ_conditions and SQ_dark_sol.
% these are properties of the set of all defects and can differ at each ij
all_Tfreezes_ij_holder = zeros(dummy_conditions.SQ_num_chardists, dummy_conditions.SQ_num_Trates, max_possible_SQ_T_values);
all_kBTfreezes_ij_holder = zeros(dummy_conditions.SQ_num_chardists, dummy_conditions.SQ_num_Trates, max_possible_SQ_T_values);

% for each ij, there will be this many temperatures we actually calculate
% at
num_unique_temps_this_ij = zeros(dummy_conditions.SQ_num_chardists, dummy_conditions.SQ_num_Trates);


%% build up the solution holder variable SQ_dark_sol
% These variables have well-defined sizes independnet of how many T's we calculate at
SQ_dark_sol.defect_names = dummy_defects.defect_names;
SQ_dark_sol.chargestate_names = dummy_defects.chargestate_names;
SQ_dark_sol.SQ_Trates = dummy_conditions.SQ_Trates;
SQ_dark_sol.SQ_chardists = dummy_conditions.SQ_chardists;

% Becasue there may be from 1 to a max number of Tfreezes, we need to make these ones here big enough to hold all possible cases
% this list is scalar values for each chardist/Trate pair so 3D arrays
SQ_dark_sol.charge_bal_err = zeros(dummy_conditions.SQ_num_chardists, dummy_conditions.SQ_num_Trates, max_possible_SQ_T_values);  %we want to save the errors for each T we do a calc at, so these holders need to be big enough
SQ_dark_sol.element_bal_err = zeros(dummy_conditions.SQ_num_chardists, dummy_conditions.SQ_num_Trates, max_possible_SQ_T_values);
SQ_dark_sol.tot_bal_err = zeros(dummy_conditions.SQ_num_chardists, dummy_conditions.SQ_num_Trates, max_possible_SQ_T_values);
SQ_dark_sol.Nd = dummy_conditions.Nd * ones(dummy_conditions.SQ_num_chardists, dummy_conditions.SQ_num_Trates, max_possible_SQ_T_values);
SQ_dark_sol.Na = dummy_conditions.Na * ones(dummy_conditions.SQ_num_chardists, dummy_conditions.SQ_num_Trates, max_possible_SQ_T_values);
SQ_dark_sol.EFn = zeros(dummy_conditions.SQ_num_chardists, dummy_conditions.SQ_num_Trates, max_possible_SQ_T_values);
SQ_dark_sol.EFp = zeros(dummy_conditions.SQ_num_chardists, dummy_conditions.SQ_num_Trates, max_possible_SQ_T_values);
SQ_dark_sol.n = zeros(dummy_conditions.SQ_num_chardists, dummy_conditions.SQ_num_Trates, max_possible_SQ_T_values);
SQ_dark_sol.p = zeros(dummy_conditions.SQ_num_chardists, dummy_conditions.SQ_num_Trates, max_possible_SQ_T_values);
if dummy_conditions.sth_flag ==1
    SQ_dark_sol.STH1 = zeros(dummy_conditions.SQ_num_chardists, dummy_conditions.SQ_num_Trates, max_possible_SQ_T_values);
    SQ_dark_sol.STH2 = zeros(dummy_conditions.SQ_num_chardists, dummy_conditions.SQ_num_Trates, max_possible_SQ_T_values);
end

% these two need to be 4D arrays - hold all the chargestates and defects at
% each computed temperature for each i,j
SQ_dark_sol.chargestates = zeros(dummy_conditions.SQ_num_chardists, dummy_conditions.SQ_num_Trates, max_possible_SQ_T_values, dummy_defects.num_chargestates);  % char. dimensions in rows, defects or chargestates in columns.
SQ_dark_sol.defects = zeros(dummy_conditions.SQ_num_chardists, dummy_conditions.SQ_num_Trates, max_possible_SQ_T_values, dummy_defects.num_defects);

%% end initializing SQ_dark_sol



% To do later: detect if any elements have a fixed number constraint - this will mean
% that the defects containing that element may end up with conflicting
% constraints when we go to freeze all of them




max_num_actual_unique_Ts_for_all_ij = 0;  % this counts the max number of T's we calculate at over all the i,j pairs so we can trim zeros off the arrays after calc is done (trim the fat)


%%%%%  Generalized Quenching calcualtion start
%% loop over the characteristic distances and cooling rates
for i = 1:working_conditions.SQ_num_chardists

    for j = 1:working_conditions.SQ_num_Trates
        disp('Please be patient, SQ calc in progress')

        % reset all defects fixed/open status to that in the original
        % conditions variable (so all are open with 0 in their fixed flag except those specified to be
        % fixed for all T in the main conditions variable).
        working_conditions.fixed_defects = dummy_conditions.fixed_defects;
        working_conditions.fixed_defects_concentrations = dummy_conditions.fixed_defects_concentrations;
        working_conditions.indices_of_fixed_defects = dummy_conditions.indices_of_fixed_defects;  %this is the starting point, we need to update this as we go
        working_conditions.num_fixed_defects = dummy_conditions.num_fixed_defects;


        % k loos over all the defects in the model
        defects_Tfreeze_this_ij = zeros(1,dummy_defects.num_defects);
        for k = 1:dummy_defects.num_defects

            % compute the Tfreezes for each defect and round them to nearest degree K

            % question for Aadi and Mike: should we jsut hard code the T interval to be
            % 1 K and thus not take SQ_T_step in as input?  Or should we modify this
            % round operation to round to the nearest SQ_T_interval?
            % also, we should round off Tstart and Tstop temperaturs to the same
            % temperature points - like we dont want someone to input 298.15 K as the
            % stop, we want to use 300 K instead.

            Do = dummy_defects.cs_Do(k) * 

            defects_Tfreeze_this_ij(k) = round(SQ_find_Tfreeze_diffusion(working_conditions.SQ_lin_or_exp, working_conditions.SQ_T_start, working_conditions.SQ_T_end, working_conditions.SQ_chardists(i), working_conditions.SQ_Trates(j), dummy_defects.cs_Do(k), dummy_defects.cs_Emigration(k)));

            % check if the calculated Tfreeze value is outside of the Tmax to Tmin range, set it to the nearby boundary value if it is.
            % additionally, if the defect is frozen for all temperatures, set its Tfreeze to the max one
            if dummy_conditions.fixed_defects(k)==1   % if the defect is always frozen.  this case comes first so will override the other possibilities in the elseifs below.  Yes it calls the original conditions variable not the local working_conditions
                Tfreeze_ij_for_defect_k_holder(i,j,k) = dummy_conditions.SQ_T_start;
            elseif defects_Tfreeze_this_ij(k) >= dummy_conditions.SQ_T_start
                Tfreeze_ij_for_defect_k_holder(i,j,k) = dummy_conditions.SQ_T_start;
            elseif defects_Tfreeze_this_ij(k) <= dummy_conditions.SQ_T_end
                Tfreeze_ij_for_defect_k_holder(i,j,k) = dummy_conditions.SQ_T_end;
            elseif (defects_Tfreeze_this_ij(k) < dummy_conditions.SQ_T_start) && (defects_Tfreeze_this_ij(k) > dummy_conditions.SQ_T_end)
                Tfreeze_ij_for_defect_k_holder(i,j,k) = defects_Tfreeze_this_ij(k);  % this is the default case where the calculated value is in the boundaries so set Tfreeze to it
            else
                error('something weird about calc_Tfreeze_ijk and the high and low T values for the SQ calc')
            end
        end  % end loop over k.  Now each defect has been assigned a Tfreeze for this ij and these are in Tfreeze_ij_for_defect_k_holder


        %To DO check if need to transpose, may be row vector already)

        % determine the list of unique freezing temperatures, agmented with the evenly-spaced list from Thigh to Tlow to help Kroger follow the
        % solution continuously rather than havign a hard time if the Tfreeze list is sparse.  Find the unique T values and sort them descending.
        % this list applies to the set of all the defects for this ij
        unique_temps_this_ij = sort(unique([defects_Tfreeze_this_ij T_reg_spaced]),'descend');   % we need to know how many values are in this to assign it to the main variable which has more spaces (so it will have a bunch of zeros after these are inserted)
        num_unique_temps_this_ij(i,j) = numel(unique_temps_this_ij);
        all_Tfreezes_ij_holder(i,j,1:num_unique_temps_this_ij(i,j)) = unique_temps_this_ij;  % all unique freezing temps including reg spaced ones.  transpose to rowvector
        all_kBTfreezes_ij_holder(i,j,1:num_unique_temps_this_ij(i,j)) = SQ_conditions.kB * all_Tfreezes_ij_holder(i,j,1:num_unique_temps_this_ij(i,j));


        % if this ij yields the largest set of T's we need to calculate at
        % so far, update the max_num_unique_Ts_for_all_ij.  After loop over
        % all ij's is doen, we can toss out any extra all-zero parts of the
        % holder variables
        if num_unique_temps_this_ij(i,j) > max_num_actual_unique_Ts_for_all_ij
            max_num_actual_unique_Ts_for_all_ij = num_unique_temps_this_ij(i,j);
        end

        % now unique_temps_this_ij AND all_Tfreezes_ij_holder(i,j,:) contain all the unique
        % freezing temperatures for this ij.  Indexing into the 3D array
        % can be a pain so keep the other one too


        %% start the actual SQ calc for this i,j pair of chardist and Trate
        last_Ef_fixed_mu_sol_vec = [];  % make these two variables but empty for now
        second_last_Ef_fixed_mu_sol_vec = [];


        for l = 1:num_unique_temps_this_ij(i,j)   % for each combination of char_dist and Trate, we need to go down the list of unique Tfreezes and do a sequential calc.

            % X_conditions.T_equilibrium, etc are the variables that the
            % main computation engine uses, so set these here one at a time
            working_conditions.T_equilibrium = unique_temps_this_ij(l);
            working_conditions.num_T_equilibrium = 1;
            working_conditions.kBT_equilibrium = working_conditions.kB * working_conditions.T_equilibrium;

            this_Tfreeze_index = find(dummy_conditions.T_equilibrium == working_conditions.T_equilibrium);  % find the index in the big look up table spaced 1 K apart for this current T we need to calc at

            working_conditions.muT_equilibrium = dummy_conditions.muT_equilibrium(this_Tfreeze_index,:); % muT_equilibrium has as many columns as elements in the model
            working_conditions.NcT_equilibrium = dummy_conditions.NcT_equilibrium(this_Tfreeze_index);
            working_conditions.NvT_equilibrium = dummy_conditions.NvT_equilibrium(this_Tfreeze_index);
            working_conditions.EgT_equilibrium = dummy_conditions.EgT_equilibrium(this_Tfreeze_index);
            working_conditions.EcT_equilibrium = dummy_conditions.EcT_equilibrium(this_Tfreeze_index);
            working_conditions.EvT_equilibrium = dummy_conditions.EvT_equilibrium(this_Tfreeze_index);

            % TO DO
            % Handle set up the ability to specify fixed defect concentrations that
            % vary vs T
            % if strcmp(dummy_conditions.T_dep_fixed_defect_flag,'On')
            %     if size(dummy_conditions.fixed_defects_concentrations,1)==Tloop_conditions.num_T_equilibrium
            %         Tloop_conditions.fixed_defects_concentrations = dummy_conditions.fixed_defects_concentrations(i1,:);          % take the i1 row of the fixed defects variable from the original conditions file and assign that row in the Tloop_conditions
            %         Tloop_conditions.fixed_defects_concentrations = dummy_conditions.fixed_defects_concentrations';  % rotate this to a column vector (doing it in 2 steps to make it obvious we did it)
            %     else
            %         error('If T_dep_fixed_defect_flag is ON, then conditions.fixed_defect_concentrations has to be a vector with an entry for each T')
            %     end
            % elseif strcmp(dummy_conditions.T_dep_fixed_defect_flag,'Off')
            %     if size(dummy_conditions.fixed_defects_concentrations,2)~=1
            %         error('If T_dep_fixed_defect_flag is OFF, then conditions.fixed_defect_concentrations has to be a single value (scalar) appropriate for all temperatures')
            %     end
            % end

            % handle T dependent dH values here
            %     if strcmp(dummy_conditions.T_dep_dH,'On')
            %         Tloop_defects.dH = dummy_defects.dHoT(i1,:);
            %
            %     elseif strcmp(dummy_conditions.T_dep_dH,'Off')
            %         %% do nothing different
            %     else
            %         error('T dependent fixed defects and or T dependent dHo flags must be on or off')
            %     end
            %


            % determine the indices of the defects that will freeze at the current Tfreeze.  This can't exceed the numebr of defects.    There may be none for Tend and Tstart, but we need to do these Temperatures anyway.
            defects_freezing_now_index = find(defects_Tfreeze_this_ij == working_conditions.T_equilibrium);
            num_defects_freezing_now = numel(defects_freezing_now_index);
            % this sometimes gives an empty 0x1 so just in case do these
            % also so that nothing weird happens
            if isempty(defects_freezing_now_index)
                % defects_freezing_now_index = []; %this is redundant with
                % the if condition
                num_defects_freezing_now = 0;
            end



            % since we are sending one temperature at a time to the main
            % calculation engine, we keep track of the last few solutions
            % and supply a guess for the Ef_mu_guess for Ef and the mu for
            % any elements with fixed concentrations.  the calc engine does
            % this automatically when T_equilibrium is a vector, but in the
            % SQ calc we are sending T values one at a time to it.  So have
            % to build in the same functionality here.

            % the main calculation at each T is done inside this if loop
            if l==1
                working_conditions.Ef_fixed_mu_guesses_supplied_flag = 'Off';
                working_conditions.Ef_fixed_mu_guesses = [];
                % do the actual calc for this T - no guesses given for
                % Ef_mu_vec
                [partial_SQ_sol] = defect_equilibrium_dark(working_conditions, dummy_defects);   % the SQ solution is built from a whole lot of equilibrium solutions with progressively more defects frozen (thus this is a partial solution only)
                [partial_fullquench_dark_sol] = defect_fullquench_dark(partial_SQ_sol, working_conditions, dummy_defects);
                last_Ef_fixed_mu_sol_vec = [partial_SQ_sol.EFn partial_SQ_sol.mu(find(working_conditions.fixed_conc_flag))];  % make a vector of the Ef_fixed_mu_vec solution for this temeperature
            elseif l==2
                working_conditions.Ef_fixed_mu_guesses_supplied_flag = 'On';
                working_conditions.Ef_fixed_mu_guesses = last_Ef_fixed_mu_sol_vec;
                % guess this time is the Ef_mu_vec that was the solution
                % for l=1
                [partial_SQ_sol] = defect_equilibrium_dark(working_conditions, dummy_defects);   % the SQ solution is built from a whole lot of equilibrium solutions with progressively more defects frozen (thus this is a partial solution only)
                [partial_fullquench_dark_sol] = defect_fullquench_dark(partial_SQ_sol, working_conditions, dummy_defects);
                second_last_Ef_fixed_mu_sol_vec = last_Ef_fixed_mu_sol_vec;
                last_Ef_fixed_mu_sol_vec = [partial_SQ_sol.EFn partial_SQ_sol.mu(find(working_conditions.fixed_conc_flag))];  % make a vector of the Ef_fixed_mu_vec solution for this temeperature
            elseif l>=3
                working_conditions.Ef_fixed_mu_guesses_supplied_flag = 'On';
                % guess this time is average of last answer and the first
                % order estimate from derivative of last two answers
                working_conditions.Ef_fixed_mu_guesses = 0.5*(last_Ef_fixed_mu_sol_vec + (second_last_Ef_fixed_mu_sol_vec - last_Ef_fixed_mu_sol_vec)/(dummy_conditions.T_equilibrium(l-2)-dummy_conditions.T_equilibrium(l-1))*(dummy_conditions.T_equilibrium(l)-dummy_conditions.T_equilibrium(l-1)));
                [partial_SQ_sol] = defect_equilibrium_dark(working_conditions, dummy_defects);   % the SQ solution is built from a whole lot of equilibrium solutions with progressively more defects frozen (thus this is a partial solution only)
                [partial_fullquench_dark_sol] = defect_fullquench_dark(partial_SQ_sol, working_conditions, dummy_defects);
                second_last_Ef_fixed_mu_sol_vec = last_Ef_fixed_mu_sol_vec;
                last_Ef_fixed_mu_sol_vec = [partial_SQ_sol.EFn partial_SQ_sol.mu(find(working_conditions.fixed_conc_flag))];  % make a vector of the Ef_fixed_mu_vec solution for this temeperature
            end
            % completed the calc at this temp, store what needs to be
            % stored

            % build up the solution structure vs ij for this value
            % of l.  4D arrays for chargestates and defects
            SQ_dark_sol.chargestates(i,j,l,:) = partial_SQ_sol.chargestates;
            SQ_dark_sol.defects(i,j,l,:) = partial_SQ_sol.defects;

            % 3D arrays
            SQ_dark_sol.EFn(i,j,l) = partial_SQ_sol.EFn;
            SQ_dark_sol.EFp(i,j,l) = partial_SQ_sol.EFp;
            SQ_dark_sol.n(i,j,l) = partial_SQ_sol.n;
            SQ_dark_sol.p(i,j,l) = partial_SQ_sol.p;
            SQ_dark_sol.Nd(i,j,l) = partial_SQ_sol.Nd;
            SQ_dark_sol.Na(i,j,l) = partial_SQ_sol.Na;
            if dummy_conditions.sth_flag == 1
                SQ_dark_sol.STH1(i,j,l) = partial_SQ_sol.STH1;
                SQ_dark_sol.STH2(i,j,l) = partial_SQ_sol.STH2;
            end
            SQ_dark_sol.charge_bal_err(i,j,l) = partial_SQ_sol.charge_bal_err;  % there will be an entry here for the number of calcs we need to do for each ij
            SQ_dark_sol.element_bal_err(i,j,l) = partial_SQ_sol.element_bal_err;
            SQ_dark_sol.tot_bal_err(i,j,l) = partial_SQ_sol.tot_bal_err;
            SQ_dark_sol.n_quench(i,j,l) = partial_fullquench_dark_sol.n;
            SQ_dark_sol.p_quench(i,j,l) = partial_fullquench_dark_sol.p;

            % now we move the defects that froze at this T over to the list of fixed
            % defects
            if num_defects_freezing_now >= 1  % detect if this is a temp when 1 or more defects froze
                for m = 1:num_defects_freezing_now %loop over the ones that froze
                    working_conditions.num_fixed_defects = working_conditions.num_fixed_defects + 1;
                    working_conditions.fixed_defects(defects_freezing_now_index(m)) = 1;   % mark each freezing defect as frozen
                    working_conditions.fixed_defects_concentrations(defects_freezing_now_index(m)) = partial_SQ_sol.defects(defects_freezing_now_index(m));   % save all their concentrations
                end  % end m loop over the defects freezing now
                working_conditions.indices_of_fixed_defects = find(working_conditions.fixed_defects==1);  % update the list of fixed defects afther these ones freeze
            elseif num_defects_freezing_now==0  % this can occur for the Tmax and Tmin but we calculate at them anyway
                disp(strcat('no defects froze at_',num2str(working_conditions.T_equilibrium),'_K'))
                working_conditions.T_equilibrium
            else
                error(strcat('number of defects freezing at_',num2str(working_conditions.T_equilibrium),'_should be 0 or greater'))
            end

        end   %end l loop over unique T freezes

        % so now we have just done the equilibrium calc at the lowest
        % (final) temperature, usually 300 K.
        % the partial_SQ_sol variable has the current concentrations of all
        % chargestates and defects in it (which includes all the frozen
        % defect concentrations)

        % at this point, all of the defects should be marked as frozen, meaning we have taken care of them,  and each
        % should have a concentration assigned.  Technically, some may not be
        % frozen at the min temp (300K) but since we have done an equilibrium
        % calc there with the others frozen, it gives the right answer at this T for all
        % the defects and chargestates.

        % why dont we leave some open at thelowest T?  Just to make it easy
        % to check if we took care of all the defects in the list and didnt
        % forget to look at any of them.

        % after looping over all of the Tfreezes, we have the full SQ solution
        % for the current characteristic dimension and Trate so we can fill in a row of
        % the output variables

        z = numel(unique_temps_this_ij);
        n = reshape(SQ_dark_sol.n(i,j,:),1,[]);
        p = reshape(SQ_dark_sol.p(i,j,:),1,[]);
        plot(unique_temps_this_ij, log10(n(1:z)),'ro')
        plot(unique_temps_this_ij, log10(n(1:z)),'bo')

        disp(strcat('Done with SQ calc for chardist=',num2str(working_conditions.SQ_chardists(i)),'and_','Trate=',working_conditions.SQ_Trates(j)))

    end  % end loop over j Trates
end % end loop over i chardists

ylim([12 19])


% % SQ_conditions.T_equilibrium =
% % SQ_conditions.kBT_equilibrium =
% % SQ_conditions.muT_equilibrium =
% % SQ_conditions.NcT_equilibrium =
% % SQ_conditions.NvT_equilibrium =
% % SQ_conditions.EgT_equilibrium =
% % SQ_conditions.EcT_equilibrium =
% % SQ_conditions.EvT_equilibrium =
% %



% trim off the fat from the holder variables in the solution (we
% preallocated these at the max possible size they could take on in terms
% of number of possible unique T's but maybe we didnt use all the space
% since for each ij there is a good chance that multiple defects will
% freeze at the same T

% % % % % % % % % % % % % % % max_num_actual_unique_Ts_for_all_ij - this one tells how many we need to keep - the extras we can trim off newxt.
% % % % % % % % % % % % % % %
% % % % % % % % % % % % % % %  num_unique_temps_this_ij - this variable tells you how many actual temperatures we calculated at for each ij.

% for each ij, how many temps were actullay used?  store the final
% num_unique_temps_this_ij matrix to help plot things later
SQ_dark_sol.l_value_for_defects_cs = num_unique_temps_this_ij;

% 4D arrays for chargestates and defects
SQ_dark_sol.chargestates = SQ_dark_sol.chargestates(:,:,1:max_num_actual_unique_Ts_for_all_ij,:);
SQ_dark_sol.defects = SQ_dark_sol.defects(:,:,1:max_num_actual_unique_Ts_for_all_ij,:);

% 3D arrays
SQ_dark_sol.EFn = SQ_dark_sol.EFn(:,:,1:max_num_actual_unique_Ts_for_all_ij);
SQ_dark_sol.EFp = SQ_dark_sol.EFp(:,:,1:max_num_actual_unique_Ts_for_all_ij);
SQ_dark_sol.n = SQ_dark_sol.n(:,:,1:max_num_actual_unique_Ts_for_all_ij);
SQ_dark_sol.p = SQ_dark_sol.p(:,:,1:max_num_actual_unique_Ts_for_all_ij);
SQ_dark_sol.Nd = SQ_dark_sol.Nd(:,:,1:max_num_actual_unique_Ts_for_all_ij);
SQ_dark_sol.Na = SQ_dark_sol.Na(:,:,1:max_num_actual_unique_Ts_for_all_ij);
SQ_dark_sol.n_quench = SQ_dark_sol.n_quench(:,:,1:max_num_actual_unique_Ts_for_all_ij);
SQ_dark_sol.p_quench = SQ_dark_sol.p_quench(:,:,1:max_num_actual_unique_Ts_for_all_ij);

if dummy_conditions.sth_flag == 1
    SQ_dark_sol.STH1 = SQ_dark_sol.STH1(i,j,1:max_num_actual_unique_Ts_for_all_ij);
    SQ_dark_sol.STH2 = SQ_dark_sol.STH2(i,j,1:max_num_actual_unique_Ts_for_all_ij);
end
SQ_dark_sol.charge_bal_err = SQ_dark_sol.charge_bal_err(:,:,1:max_num_actual_unique_Ts_for_all_ij);
SQ_dark_sol.element_bal_err = SQ_dark_sol.element_bal_err(:,:,1:max_num_actual_unique_Ts_for_all_ij);
SQ_dark_sol.tot_bal_err = SQ_dark_sol.tot_bal_err(:,:,1:max_num_actual_unique_Ts_for_all_ij);

SQ_dark_sol.Tfreeze_ij_for_defect_k_holder = Tfreeze_ij_for_defect_k_holder;
% % % % % % % % % % % % % SQ_dark_sol.T_equilibrium =
% % % % % % % % % % % % % SQ_dark_sol.num_T_equilibrium =


% % % % % % % % % % % fill in the SQ_conditions for a complete record of the calculation, and fill in te SQ_darksol



end   %%%%%%%%%  end of the actual SQ routine
%% function that computes the Tfreeze for each defect based on diffusion length for cooling at a constant rate from Tmax to Tfreeze
% T(t) = Tmax - rate*time

% Arrhenius equation asssumed for D(T) Do.*exp(-Ea./kB.*T(t)).
% expression originally derived in mathematica.
% The definition of expint in Matlab is % different from Mathematica's ExpIntEi.
%  Expint(z) in matlab = integral(1./t.*exp(-t) from z to infinity.
%  Mathematica has the lower limit of integtal as -z while matlab has +z.
% So mathematica expint(x) = -expint(-x) in matlab.

% this function will be called one defect at a time.  
% variable types for inputs: 
% lin/exp = string
% current_chardist = scalar
% current_Trates = scalar
% Tmax, Tmin = scalar
% Do, Ea = scalar
% Output T_freeze = vector



function [T_freeze] = SQ_find_Tfreeze_reaction(lin_or_exp, Tmax, Tmin, current_chardist, current_Trate, Do, Ea)

kB_eV = 8.617333262e-5;  % eV/K
% kB_J = 1.380649e-23; % J/K
% q = 1.60217663e-19;  % J/eV

Tfreeze_int = (Tmax-Tmin)/100;  % this anticipates Tmax and Tmin might be numbers like 1800 K and 300 K.  So we generate a grid of 100 K increments 
Tfreeze_list = Tmax:-Tfreeze_int:Tmin;   % coarse list of freezing temperatures across the whole temperature rnage.  essentially guesses on Tfreeze to bracket the actual one.   note the order from high T to low T 

## what is the criterion here (not a distance but a rxn extent or fractional change in number or something like that)
xdiff_list = DiffLength(Tfreeze_list);  
## 


if ~isreal(xdiff_list)
    error('Some values of xdiff_guesses were complex, so negative valies for difflength^2 generated')
end


if current_chardist ==0  % this is the case for the first char_dist in the list.  This has to come before the other cases so it catches it
    T_freeze = Tmin;

elseif xdiff_list(1) <= current_chardist    % case where diffusion length is never long enough, or just long enough at Tmax - set Tf to Tmax
    T_freeze = Tmax;

elseif xdiff_list(end) >= current_chardist   % case where diffusion length is alwasy long enough - set Tf to Tmin
    T_freeze = Tmin;

elseif (xdiff_list(1) > current_chardist) && (xdiff_list(end) < current_chardist)  % normal case where Tf falls between Tmax and Tmin.  Note Tfreeze_list is decreasing from 1:end
    [T_guess] = T_freeze_guess(Tfreeze_list, xdiff_list, current_chardist);   % this line figures out two starting temperatures for fzero to start at.  
    [T_freeze,~,exitflag,~] = fzero(@err_func,T_guess);  %tis finds the actual Tfreeze
    if exitflag<0
        exitflag
        error('fzero terminated with unusual conditions while trying to find tfreeze-time to debug')
    end

elseif current_chardist < 0
    current_chardist
    xdiff_list
    error('current_chardist is negative')
else
    current_chardist
    xdiff_list
    error('something strange about xdiff - examine it')
end





%%%%%%%%%%%%%%% nested subroutines %%%%%%%%%%%%%%%%%%%%%%%%%%%

%% function that guesses starting points for fsolve
    function [T_guess] = T_freeze_guess(Tfreeze_list, xdiff_list, char_dist)

        edge_index = find(diff(sign(xdiff_list-char_dist))~=0);   %  this finds the rising/falling edge where (xdiff - chardist) changes sign.  Note the diff operation shortens the vector by one entry.  

        if numel(edge_index)==1   % we get 1 entry if the sign changes but we dont have the exact Tfreeze in our list - normal case
            T_guess = [Tfreeze_list(edge_index) Tfreeze_list(edge_index+1)];   %create a 2 entry guess

        elseif numel(edge_index)==2 && diff(edge_index)==1  %this can happen if the Tfreeze is actually coincidentally in our list so produces a zero in the diff - it will be the 2nd value if this happens.
            T_guess = Tfreeze_list(edge_index);  %take both entries 

        elseif numel(edge_index)==2 && diff(edge_index)~=1 || numel(edge_index)>2  % an edge case could occur if somehow difflength - char_dist is non-monotonic and has 2 or more solutions
            xdiff_list
            char_dist
            error('difflength - char dist seems like it may not be monotonic')
        else
            xdiff_list
            char_dist
            error('something strange about the function vs Tmin - solutions may not be valid')
        end

    end  %%%% end T_freeze_guess


    function [xdiff_list] = DiffLength(Tfreeze_list)    % function for the equation to give difflength as function of the vector of starting guess Tfreeze values
        % in this function, the arguments of special functions can be in eV, but
        % the other factors all need to be watched for units
        % Trates is a vector, T_freeze_guess is scalar, Do and Ea are
        % scalars
        % Notes: Herein we assume z is a real number.  ExpintEi in mathematica is (-1)*integral of (e-t)/t dt
        % from -z to infinity. Define this as Ei_Mathematica(z).  To convert to expint() function in matlab:   Ei_Mathematica = -expint(-z).  
        % Also there are differences in how the two programs define the
        % incomplete Gamma function.  Gamma_Mathematica[0,z] =
        % -Ei_Mathematica(-z) which turns into expint(z) in Matlab. 


        % the solution for linear cooling in Mathematica follows
        % precompute some quantities
        Xa = Ea/kB_eV;
        Xf = Xa./Tfreeze_list;
        Xm = Xa/Tmin;
        if strcmp(lin_or_exp,'lin')
            xdiff_list = sqrt( Do/current_Trate*( Tfreeze_list.*exp(-Xf) - Tmin*exp(-Xm) - Xa*real(expint(Xf)) + Xa*real(expint(Xm))) );  %triple checked these formulae from Mathematica MAS 1/27/26
        elseif strcmp(lin_or_exp,'exp')
            xdiff_list = sqrt( Do/current_Trate * Tfreeze_list * real(expint(Xf)) );  % checked by MS 1-26-26.  This is an approximate solution where T goes all the way to 0 Kelvin but not much will happen below room temp usually so it's pretty OK.  See Mathematica notebook from 1/26/2026 for details
        else
            error('exp_or_lin must be a string saying exp or lin')
        end



    end




    function [err] = err_func(Tfreeze)    % function to give the difflength - char_dist.  This will send in a scalar T_freeze
        err = DiffLength(Tfreeze) - current_chardist;
    end



    function [Keq] = find_Keq(defects, rxns, kBT)
        Keq = exp(-rxns.delta_G/kBT);
    end




    function [kf, kb] = rate_constants(rxns, defects, kBT)
        % generate the current rate constants for all reactions given
        % temperature and defect concentraitons (accounting for
        % defetc-mediated diffusion).  Model is diffusion-limited capture
        % model.  the rxns variable can hold different values for each
        % complex - this will allow for example the inclusion of multiple reactions
        % forming "excited states" of complexes - meaning the same two
        % defects separated by different distances or on different pairs of
        % sites and thus having different binding energies realtive to the separated defects.  In that case the capture radius should be the
        % inter-defect spacing itself in that config. If configs of the
        % complex are lumped together, then the cpature radius should be
        % computed from the idea of a capture event horizon: kBT = |coulomb potential(r)| - the radius below
        % which the formation of the bound pair in ground state as next-door neighbors is inevitable
        cs_D = diffusion_consts(defects, current_concentrtions, kBT);  % calc the D values for all defects at this T
        Keq = find_Keq(defects, rxns, kBT);
        for i=1:rxns.num_rxns
            kf(i) = 4*pi*rxns.capture_radius(i)*(cs_D(rxns.reactants(i,1)) + cs_D(rxns.reactants(i,2)))*exp(-rxns.E_activation/kBT); % calculate kf including possibility of an activation barrier for the reaction
        end
        kb = kf./Keq;  % element-wide vector divide, compute 
    end


        function [D] = diffusion_consts(defects, current_concentrations, kBT)
            % take in the temperature, defect database, and current
            % concentraitons of defect and spit out the diffusion const for
            % each chargestate
            D = defects.

        end
       




    function dC/dt  = rates(chargestate conc, n, p, sth1, sth2, reactions)





        reactions = database of reactions amongst a set of chargestates
        - have to list all possible reactions amongst chargestates and carriers (in general)
        - only binary reactions
        - list reactants and products in sorted order with no duplicates (check this upon forming the rxn database)
        - will be a Nx2 matrix
        - check that no reactans show up as products in same reaction
        reactions are mass action kinetics with exponents same as stoich for the reaction











end   %%% end main function
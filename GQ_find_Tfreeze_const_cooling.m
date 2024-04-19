%% function that computes the Tfreeze for each defect based on diffusion length for cooling at a constant rate from Tmax to Tfreeze
% T(t) = Tmax - rate*time

% Arrhenius equation asssumed for D(T) Do.*exp(-Ea./kB.*T(t)).
% expression originally derived in mathematica.
% The definition of expint in Matlab is % different from Mathematica's ExpIntEi.
%  Expint(z) in matlab = integral(1./t.*exp(-t) from z to infinity.
%  Mathematica has the lower limit of integtal as -z while matlab has +z.
% So mathematica expint(x) = -expint(-x) in matlab.

function [T_freeze] = GQ_find_Tfreeze_const_cooling(lin_or_exp, char_dist, Tmax, Tmin, Trate, Do, Ea)

kB_eV = 8.617333262e-5;  % eV/K
% kB_J = 1.380649e-23; % J/K
% q = 1.60217663e-19;  % J/eV

Tfreeze_int = -(Tmax-Tmin)/200;
Tfreeze_list = Tmax:Tfreeze_int:Tmin;   %note the order from high T to low T
xdiff_list = DiffLength(Tfreeze_list);

if ~isreal(xdiff_list)
    error('Some values of xdiff were complex, so negative valies for difflength^2 generated')
end

% figure(3)
% clf
% hold on
% plot(Tfreeze_grid,xdiff)
% plot(Tfreeze_grid,char_dist*ones(size(Tfreeze_grid)))
% %         ylim([-kB_eV kB_eV])


if char_dist==0  % this is the case for the first char_dist in the list.  This has to come before the other cases so it catches it
    T_freeze = Tmin;

elseif xdiff_list(1)<=char_dist    % case where diffusion length is never long enough, or just long enough at max value - set Tf to Tmax
    T_freeze = Tmax;

elseif xdiff_list(end)>=char_dist   % case where diffusion length is alwasy long enough - set to Tmin
    T_freeze = Tmin;

elseif xdiff_list(1)>char_dist && xdiff_list(end)<char_dist  % normal case where Tf falls between Tmax and Tmin
    [T_guess] = T_freeze_guess(Tfreeze_list, xdiff_list, char_dist);
    [T_freeze,~,exitflag,~] = fzero(@err_func,T_guess);
    if exitflag<0
        exitflag
        error('fzero terminated with unusual conditions while trying to find tfreeze-time to debug')
    end

elseif char_dist<0
    char_dist
    xdiff_list
    error('char_dist is negative')
else
    char_dist
    xdiff_list
    error('something strange about xdiff - examine it')
end





%%%%%%%%%%%%%%% nested subroutines %%%%%%%%%%%%%%%%%%%%%%%%%%%

%% function that guesses starting points for fsolve
    function [T_guess] = T_freeze_guess(Tfreeze_list, xdiff_list,char_dist)

        edge_index = find(diff(sign(xdiff_list-char_dist))~=0);   %  this finds the rising/falling edge where (xdiff - chardist) changes sign.

        if numel(edge_index)==1   % we get 1 entry if the sign changes but we dont have the exact Tfreeze in our list - normal case
            T_guess = Tfreeze_list(edge_index);

        elseif numel(edge_index)==2 && diff(edge_index)==1  %this can happen if the Tfreeze is actually coincidentally in our list - it will be the 2nd value if this happens.
            T_guess = Tfreeze_list(edge_index(2));

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




    function [xdiff_list] = DiffLength(Tfreeze_list)    % function for the equation to give difflength as function of trial Tfreeze values
        % in this function, the arguments of special functions can be in eV, but
        % the other factors all need to be watched for units
        if strcmp(lin_or_exp,'lin')
            xdiff_list = sqrt(  Do/Trate*(Tfreeze_list.*exp(-Ea./(kB_eV*Tfreeze_list)) - Tmin*exp(-Ea/(kB_eV*Tmin)) + Ea/kB_eV*(real(expint(Ea/(kB_eV*Tmin))) - real(expint(Ea./(kB_eV*Tfreeze_list)))) ));   %expression checked MAS & AA 4/11/23 - updated to keep rate constant (not to) during integration
        elseif strcmp(lin_or_exp,'exp')
            xdiff_list = sqrt(Do*Tmax/Trate*real(expint(Ea./(kB_eV*Tfreeze_list))));  %expression checked MAS & AA 4/11/23
        else
            error('exp_or_lin must be a string saying exp or lin')
        end
    end



    function [err] = err_func(Tfreeze)    % function to give the difflength - char_dist
        err = DiffLength(Tfreeze) - char_dist;
    end


end   %%% end main function
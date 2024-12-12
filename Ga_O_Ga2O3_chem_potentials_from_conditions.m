%% Function that returns the Ga and O chemical potentials in equilibrium with Ga2O3 for a vector of temperatures
%% 6/25/2024 - Aadi merged prior subroutines into one file so that new materials can be made using only one file.
%% In the Ga-O system, we have Ga, GaO(vap), Ga2O (s, vap), Ga2O3(s,l), and O2(g).  The vapor formed over Ga2O3 is governed by Ga2O3 = Ga2O + O2.
%% pure liquid Ga + Ga2O3 results in Ga2O evolving according to  4 Ga + 2 Ga2O3 = 4 Ga2O + O2.

function [mu_Ga, mu_O] = Ga_O_Ga2O3_chem_potentials_from_conditions(conditions)

disp('Assumption is that the specified chemical potentials are continuously enforced at constant P (no products build up or reactants deplete). This is opposed to setting these initial conditions in a closed box and letting them react.  So user make sure specified conditions are compatible and undrstood.')

if strcmp(conditions.T_dep_matrix_mu_flag,'Off')   %T independent mu values
    disp('using T-independent chemical potentials')

    if strcmp(conditions.mu_conditions_flag,'O-rich')
        mu_O = 0;
        mu_Ga = conditions.mu_Ga2O3/2;

    elseif strcmp(conditions.mu_conditions_flag,'Ga-rich')   %% for now main one should be O-rich. Next the immediate next phase will be in place of "Ga"-rich
        mu_Ga = 0;
        mu_O = conditions.mu_Ga2O3/3;
    else
        error('using T-independent mu values, conditions.mu_conditions_flag should be O-rich or Ga-rich')
    end

    % user needs to commit these to the conditions variable in the main script
    mu_Ga = ones(size(conditions.T_equilibrium,2),1)*mu_Ga;
    mu_O = ones(size(conditions.T_equilibrium,2),1)*mu_O;


elseif strcmp(conditions.T_dep_matrix_mu_flag,'On')   %% T dependent chemical potentials for the O and Ga
    disp('using T-dependent chemical potentials')

    if strcmp(conditions.mu_conditions_flag,'O2-1atm')   % equilibrium with O2 at pO2=ptot=1 atm
        conditions.pGa_vap = 0;
        conditions.pGaO_vap = 0;
        conditions.pGa2O_vap = 0;
        conditions.pO2 = 1;
        conditions.P_tot = 1;   % this overrides the total pressure if it is not 1
        conditions.P_units = 'atm';  % this scenario overrides any earlier setting of pressure units
        [G0_Ga2O3] = G0_Ga2O3_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units);   % in these ones the pressure does nothing except sizing the output array (it is still needed though)
        [G0_O2] = G0_O2_gv(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units);
        mu_O = G0_O2/2;
        mu_Ga = (G0_Ga2O3 - 3*mu_O)/2;

    elseif strcmp(conditions.mu_conditions_flag,'air-1atm')
        conditions.pGa_vap = 0;
        conditions.pGaO_vap = 0;
        conditions.pGa2O_vap = 0;
        conditions.pO2 = 0.2;
        conditions.pN2 = 0.8;
        conditions.P_tot = 1;   % this overrides the total pressure if it is not 1
        conditions.P_units = 'atm';  % this scenario overrides any earlier setting of pressure units.
        [G0_Ga2O3] = G0_Ga2O3_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units);   % in these ones the pressure does nothing except sizing the output array (it is still needed though)
        X_O2 = conditions.pO2/conditions.P_tot;
        X_N2 = conditions.pN2/conditions.P_tot;
        [G0_O2] = G0_O2_gv(conditions.T_equilibrium, conditions.P_tot, X_O2, conditions.P_units);
        mu_O = G0_O2/2;
        mu_Ga = (G0_Ga2O3 - 3*mu_O)/2;

    elseif strcmp(conditions.mu_conditions_flag,'pO2-variable')   % equilibrium with pO2 set to arbitrary value.  you must also set pO2.
        conditions.pGa_vap = 0;
        conditions.pGaO_vap = 0;
        conditions.pGa2O_vap = 0;
        if isempty(conditions.pO2)
            error('if using conditions.mu_conditions_flag = pO2-variable, user must specify a pO2')
        elseif isempty(conditions.P_units)
            error('if using conditions.mu_conditions_flag = pO2-variable, user must specify P units')
        end
        [G0_Ga2O3] = G0_Ga2O3_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units);
        X_O2 = conditions.pO2/conditions.P_tot;
        [G0_O2] = G0_O2_gv(conditions.T_equilibrium, conditions.P_tot, X_O2, conditions.P_units);  % this allows for variable pO2
        mu_O = G0_O2/2;
        mu_Ga = (G0_Ga2O3 - 3*mu_O)/2;

    elseif strcmp(conditions.mu_conditions_flag,'pGa_vap-variable')   % equilibrium with Ga_vapor.  Must also specify a pGa_vap
        conditions.pO2 = 0;
        if isempty(conditions.pGa_vap)
            error('if using conditions.mu_conditions_flag = pGa_vap-variable, user must specify a pGa_vap')
        elseif isempty(conditions.P_units)
            error('if using conditions.mu_conditions_flag =  pGa_vap-variable, user must specify P units')
        end
        conditions.pGaO_vap = 0;
        conditions.pGa2O_vap = 0;
        [G0_Ga2O3] = G0_Ga2O3_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units);
        X_Ga_vap = conditions.pGa_vap/conditions.P_tot;
        mu_Ga = G0_Ga_vap(conditions.T_equilibrium, conditions.P_tot, X_Ga_vap, conditions.P_units);
        mu_O = (G0_Ga2O3 - 2*mu_Ga)/3;

    % elseif strcmp(conditions.mu_conditions_flag,'Ga_equilib_vap_pressure')   % equilibrium vapor pressure of Ga ovr Ga
    %     conditions.pGa2O_vap = 0;
    %     conditions.pGaO_vap = 0;
    %     conditions.pO2 = 0;
    %     conditions.pGa_vap = Ga_equilib_vap_pres(conditions.T_equilibrium, conditions.P_units);
    %     if isempty(conditions.P_units)
    %         error('if using conditions.mu_conditions_flag =  Ga_equilib_vap_pressure, user must specify P units')
    %     end

    elseif strcmp(conditions.mu_conditions_flag,'Ga_ls-rich')   % equilibrium of Ga2O3 in direct contact with pure liquid or solid Ga
        conditions.pGa_vap = 0;
        conditions.pGa2O_vap = 0;
        conditions.pGaO_vap = 0;
        conditions.pO2 = 0;
        [G0_Ga2O3] = G0_Ga2O3_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units);   % pure condensed phases have X_i=1
        mu_Ga = G0_Ga_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units);
        mu_O = (G0_Ga2O3 - 2*mu_Ga)/3;

    % elseif strcmp(conditions.mu_conditions_flag,'pGa2O_vap-variable')   % equilibrium of Ga2O3 in contact with a specified p_Ga2O
    %     conditions.pGa_vap= 0;
    %     conditions.pGaO_vap = 0;
    %     conditions.pO2 = 0;
    %     if isempty(conditions.pGa2O_vap)
    %         error('if using conditions.mu_conditions_flag = pGa2O_vap-variable, user must specify a pGa2O')
    %     elseif isempty(conditions.P_units)
    %         error('if using conditions.mu_conditions_flag =  pGa2O_vap-variable, user must specify P units')
    %     end

    % elseif strcmp(conditions.mu_conditions_flag,'pGaO_vap-variable')   % equilibrium of Ga2O3 in contact with a specified p_GaO
    %     conditions.pGa_vap= 0;
    %     conditions.pGa2O_vap = 0;
    %     conditions.pO2 = 0;
    %     if isempty(conditions.pGaO_vap)
    %         error('if using conditions.mu_conditions_flag = pGaO_vap-variable, user must specify a pGaO')
    %     elseif isempty(conditions.P_units)
    %         error('if using conditions.mu_conditions_flag =  pGaO_vap-variable, user must specify P units')
    %     end

        %
        % elseif strcmp(conditions.mu_conditions_flag,'Ga2O-O2-Ga2O3_equilibrium')   % Ga2O3 in equilibrium with Ga2O(vap) and O2 vapor from decomposition of some other Ga2O3.  This is a fixed point (muO,muGa) for each T so doesnt use any pressures
        %     conditions.pGa_vap = 0;
        %     conditions.pGa2O = 0;
        %     conditions.pGaO = 0;
        %     conditions.pO2 = 0;
        %
        % elseif strcmp(conditions.mu_conditions_flag,'GaO-O2-Ga2O3_equilibrium')   % Ga2O3 in equilibrium with GaO(vap) and O2 vapor from decomposition of some other Ga2O3.  This is a fixed point (muO,muGa) for each T so doesnt use any pressures
        %     conditions.pGa_vap = 0;
        %     conditions.pGa2O = 0;
        %     conditions.pGaO = 0;
        %     conditions.pO2 = 0;

        %     elseif manual case to allow growth or etching

    else
        error('using T dependent chemical potentials, conditions.mu_conditions_flag must be one of the coded scenarios  ')
    end  % end main if elseif else structure




    %% these two lines make sure that all values are scalars, and that the total pressure is >= to the sum of reactive (not inert gas) species
    if ~isscalar(conditions.P_tot) || ~isscalar(conditions.pGa_vap) || ~isscalar(conditions.pO2) || ~isscalar(conditions.pGaO_vap) || ~isscalar(conditions.pGa2O_vap)     % check that each pressure is a scalar
        error('all pressures must scalar values  - make a loop over pressure in the script calling this if multiple P values are desired')
    end

    if conditions.P_tot<(conditions.pGa_vap + conditions.pO2 + conditions.pGaO_vap + conditions.pGa2O_vap)   % check that the total pressure (whcih could contain inert gases) is greater than or equal to sum of the specified species
        error('the total pressure must be >= sum of reactive (non inert) species pressures')
    end



end

end   % end main function






%     % work on this one - right now it doesnt actually compute the
%     % mu_Ga only the Ga vapor pressure
% elseif strcmp(conditions.mu_conditions_flag,'Ga_equilib_vap_pressure')
%     [G0_Ga2O3] = G0_Ga2O3_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units);
%     conditions.pGa_vap = Ga_equilib_vap_pres(conditions.T_equilibrium, conditions.P_units);
%     [G0_Ga_vapor] = G0_Ga_vap(conditions.T_equilibrium, conditions.P_tot, X_Ga_vap, conditions.P_units);
%
%
%     mu_Ga = conditions.pGa_vap;
%     mu_O = (G0_Ga2O3 - 2*mu_Ga)/3;
%
% elseif strcmp(conditions.mu_conditions_flag,'pGa2O_vap-variable')
%     [G0_Ga2O3] = G0_Ga2O3_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units);
%     X_Ga2O = conditions.pGa2O_vap/conditions.P_tot;
%     [G0_Ga2O] = G0_Ga2O_gv(conditions.T_equilibrium, conditions.P_tot, X_Ga2O, conditions.P_units);
%     mu_Ga = (3*G0_Ga2O - G0_Ga2O3)/4;
%     mu_O = (G0_Ga2O3 - G0_Ga2O)/2;
%
% elseif strcmp(conditions.mu_conditions_flag,'pGaO_vap-variable')
%     [G0_Ga2O3] = G0_Ga2O3_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units);
%     X_GaO = conditions.pGaO_vap/conditions.P_tot;
%     [G0_GaO] = G0_GaO_gv(conditions.T_equilibrium, conditions.P_tot, X_GaO, conditions.P_units);
%     mu_Ga = (3*G0_GaO - G0_Ga2O3);
%     mu_O = (G0_Ga2O3 - 2*G0_GaO);

%% these ones need some more thought and careful treatment becasue we need to be careful to use the stable phase of Ga2O or GaO for each T, and these could be s, l, or vap depending on T and P
%% this could also be approached by calculating for the reaction of 4Ga + 2Ga2O3 = 4Ga2O + O2 instead.
% % elseif strcmp(conditions.mu_conditions_flag,'Ga2O_vap-O2_g-Ga2O3_s-3phase')
% %     [G0_Ga2O3] = G0_Ga2O3_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units);
% %     [G0_Ga2O] = G0_Ga2O_gv(conditions.T_equilibrium, conditions.P_tot, ??, conditions.P_units);
% %     mu_Ga = 2*G0_Ga2O3 - 3*G0_Ga2O;
% %     mu_O = 2*G0_Ga2O - G0_Ga2O3;
% %
% % elseif strcmp(conditions.mu_conditions_flag,'GaO_vap-O2_g-Ga2O3_s-3phase')
% %     [G0_Ga2O3] = G0_Ga2O3_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units);
% %     [G0_GaO] = G0_GaO_gv(conditions.T_equilibrium, conditions.P_tot, ??, conditions.P_units);
% %     mu_Ga = (3*G0_GaO - G0_Ga2O3);
% %     mu_O = (G0_Ga2O3 - 2*G0_GaO);




%% subroutine for calculating vapor pressure of Ga over Ga(l)

% To do on this one: make it return both the pressure and mu values
%
% function [Ga_equilib_vap_pres] = Ga_equilib_vap_pres(T, P_units)
%     %% Gives the equilibrium vapor pressure of Ga over elemental Ga vs T
%     % T is a vector
%
%     % constants
%     kB_eV = 8.617333262e-5;
%
%     if strcmp(P_units,'atm')
%         P_ref = 1;
%     elseif strcmp(P_units,'Torr')
%         P_ref = 760;
%     elseif strcmp(P_units,'Bar')
%         P_ref = 1;
%     elseif strcmp(P_units,'Pa')
%         P_ref = 1e5;
%     else
%         error('Units of pressure must be atm, Torr, Pa, or Bar')
%     end
%
%     T = T(:);   %make T a col vactor
%
%     % compute the pressure-independent difference in G0 between the vapor and
%     % condensed phases.
%     [G0_Ga_ls] = G0_Ga_ls(T, P_ref, 1, P_units);
%     [G0_Ga_gv] = G0_Ga_gv(T, P_ref, 1, P_units);
%    Ga_equilib_vap_pres = P_ref*exp(-(G0_Ga_gv-G0_Ga_ls)./(kB_eV*T));
%
% mu_Ga =
% end
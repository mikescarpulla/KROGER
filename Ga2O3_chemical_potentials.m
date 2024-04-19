function [mu_Ga, mu_O] = Ga2O3_chemical_potentials(conditions) 

% function [mu_Ga, mu_O] = Ga2O3_chemical_potentials(conditions.T_equilibrium, conditions.P_tot, conditions.pGa_vap, P_O2, conditions.pGa2O_vap, conditions.pGaO_vap, conditions.P_units, conditions.mu_conditions_flag) 
%% Gives the chemical potentials of Ga and O in three condition
%% conditions.T_equilibrium is a vector list
%% P is a fixed value
%% conditions.P_units = bar, atm, etc
%% rich_flag = O-rich or Ga-rich or Ga2O-rich
%% O-rich solves the value for unknown mu_Ga
%% Ga-rich solves the value for unknown mu_O
%% Ga2O solves the value for unknown mu_Ga and mu_O when Ga2O3 is equilibrium with Ga2O
%% 0 < conditions.T_equilibrium < 3700

% constants
q = 1.602176634e-19;
avo = 6.0221409e+23;
kB_eV = 8.617333262e-5;

disp('Assumption is that the specified chemical potentials are continuously enforced at constant P (no products build up or reactants deplete). This is opposed to setting these initial conditions in a closed box and letting them react.  So user make sure specified conditions are compatible and undrstood.')

if strcmp(conditions.P_units,'atm')
    P_ref = 1;
elseif strcmp(conditions.P_units,'Torr')
    P_ref = 760;   
elseif strcmp(conditions.P_units,'Bar')
    P_ref = 1;
elseif strcmp(conditions.P_units,'Pa')
    P_ref = 1e5;
else
    error('Units of pressure must be atm, Torr, Pa, or Bar')
end


%%% the main calculation routines can not handle T-P 2D arrays, only
%%% vectorized for T so have to do each P separately in another loop
%%% somehow 
if isscalar(conditions.P_tot) && isscalar(conditions.pGa_vap) && isscalar(conditions.pO2) && isscalar(conditions.pGa2O_vap) && isscalar(conditions.pGaO_vap)     % check that each pressure is a scalar

    if conditions.P_tot>=(conditions.pGa_vap + conditions.pO2 + conditions.pGa2O_vap + conditions.pGaO_vap)   % check that the total pressure (whcih could contain inert gases) is greater than or equal to sum of the specified species

        if strcmp(conditions.mu_conditions_flag,'O2-1atm')
            [G0_Ga2O3] = G0_Ga2O3_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units);   % in these ones the pressure does nothing except sizing the output array (it is still needed though)
            [G0_O2] = G0_O2_gv(conditions.T_equilibrium, 1, 1, 'atm');
            mu_O = G0_O2/2;
            mu_Ga = (G0_Ga2O3 - 3*mu_O)/2;

        elseif strcmp(conditions.mu_conditions_flag,'pO2-variable')
            [G0_Ga2O3] = G0_Ga2O3_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units);   
            X_O2 = conditions.pO2/conditions.P_tot;
            [G0_O2] = G0_O2_gv(conditions.T_equilibrium, conditions.P_tot, X_O2, conditions.P_units);  % this allows for variable pO2
            mu_O = G0_O2/2;
            mu_Ga = (G0_Ga2O3 - 3*mu_O)/2;

        elseif strcmp(conditions.mu_conditions_flag,'Ga_ls-rich')
            [G0_Ga2O3] = G0_Ga2O3_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units);   % pure condensed phases have X_i=1
            [G0_Ga_ls] = G0_Ga_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units);
            mu_Ga = G0_Ga_ls;
            mu_O = (G0_Ga2O3 - 2*mu_Ga)/3;

        elseif strcmp(conditions.mu_conditions_flag,'pGa_vap-variable')
            [G0_Ga2O3] = G0_Ga2O3_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units);   
            X_Ga_vap = conditions.pGa_vap/conditions.P_tot;
            [G0_Ga_vapor] = G0_Ga_vap(conditions.T_equilibrium, conditions.P_tot, X_Ga_vap, conditions.P_units);
            mu_Ga = G0_Ga_vapor;
            mu_O = (G0_Ga2O3 - 2*mu_Ga)/3;

        elseif strcmp(conditions.mu_conditions_flag,'Ga_equilib_vap_pressure')
            [G0_Ga2O3] = G0_Ga2O3_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units);   
            conditions.pGa_vap = Ga_equilib_vap_pres(conditions.T_equilibrium, conditions.P_units);  
            mu_Ga = conditions.pGa_vap;
            mu_O = (G0_Ga2O3 - 2*mu_Ga)/3;

        elseif strcmp(conditions.mu_conditions_flag,'pGa2O_vap-variable')
            [G0_Ga2O3] = G0_Ga2O3_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units);  
            X_Ga2O = conditions.pGa2O_vap/conditions.P_tot;
            [G0_Ga2O] = G0_Ga2O_gv(conditions.T_equilibrium, conditions.P_tot, X_Ga2O, conditions.P_units);  
            mu_Ga = (3*G0_Ga2O - G0_Ga2O3)/4; 
            mu_O = (G0_Ga2O3 - G0_Ga2O)/2;

        elseif strcmp(conditions.mu_conditions_flag,'pGaO_vap-variable')
            [G0_Ga2O3] = G0_Ga2O3_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units);  
            X_GaO = conditions.pGaO_vap/conditions.P_tot;
            [G0_GaO] = G0_GaO_gv(conditions.T_equilibrium, conditions.P_tot, X_GaO, conditions.P_units);  
            mu_Ga = (3*G0_GaO - G0_Ga2O3); 
            mu_O = (G0_Ga2O3 - 2*G0_GaO);


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


% conditions.mu_conditions_flag = 'Ga2O-O2-Ga2O3-3phase';    % these two conditions specify unique mu_O and mu_Ga because they are invariant points
    % conditions.mu_conditions_flag = 'GaO-O2-Ga2O3-3phase';  

            
        else
            error('rich_flag must be O2_1atm, pO2-variable, Ga_ls_rich, pGa_vap_variable, pGa2O_vap-variable, or pGaO_vap-variable for now')
        end

    else
        error('the total pressure must be >= sum of reactive species pressures')
    end

else
    error('all pressures must be one scalar value')
end




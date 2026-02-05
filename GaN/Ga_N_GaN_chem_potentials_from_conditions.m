%% Function that returns the Ga and N chemical potentials in equilibrium with GaN for a vector of temperatures

function [mu_Ga, mu_N] = Ga_N_GaN_chem_potentials_from_conditions(conditions)

disp('Assumption is that the specified chemical potentials are continuously enforced at constant P (no products build up or reactants deplete). This is opposed to setting these initial conditions in a closed box and letting them react.  So user make sure specified conditions are compatible and undrstood.')

if strcmp(conditions.T_dep_matrix_mu_flag,'Off')   %T independent mu values
    disp('using T-independent chemical potentials')

    if strcmp(conditions.mu_conditions_flag,'N-rich')
        mu_N = 0;
        mu_Ga = conditions.mu_GaN;
    elseif strcmp(conditions.mu_conditions_flag,'Ga-rich')   %% for now main one should be O-rich. Next the immediate next phase will be in place of "Ga"-rich
        mu_Ga = 0;
        mu_N = conditions.mu_GaN;
    else
        error('using T-independent mu values, conditions.mu_conditions_flag should be Ga-rich or N-rich')
    end

    % user needs to commit these to the conditions variable in the main script
    mu_Ga = ones(size(conditions.T_equilibrium,2),1)*mu_Ga;
    mu_N = ones(size(conditions.T_equilibrium,2),1)*mu_N;

    conditions.pGa_vap = 0;



elseif strcmp(conditions.T_dep_matrix_mu_flag,'On')   %% T dependent chemical potentials from thermochemistry
    disp('using T-dependent chemical potentials')

    if strcmp(conditions.mu_conditions_flag, 'GaN_touching_Ga(l,s)')
        mu_Ga = G0_Ga_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units);
        mu_N = G0_GaN_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units) - mu_Ga;
    
    elseif strcmp(conditions.mu_conditions_flag, 'Ga_equilibrium_vapor_pressure_over_Ga_at_T')
        G_Ga_condensed = G0_Ga_ls(conditions.T_Ga, conditions.P_tot, conditions.X, conditions.P_units);
        G_Ga_vap = G0_Ga_gv(conditions.T_Ga, conditions.P_tot, 1, conditions.P_units);
        mu_Ga = conditions.P_ref*exp(-(G_Ga_vap-G_Ga_condensed)/conditions.kBT_equilibrium);
        conditions.pGa_vap = conditions.P_tot * exp(-mu_Ga/conditions.kBT_equilbrium);
        mu_N = G0_GaN_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units) - mu_Ga;

    elseif strcmp(conditions.mu_conditions_flag, 'pN2-variable') %equilibrium with pN2 set to arbitrary value. You must also set pN2.
        if isempty(conditions.p_N2)
            error('if using conditions.mu_conditions_flag = pN2-variable, user must specify a pN2')
        elseif isempty(conditions.P_units)
            error('if using conditions.mu_conditions_flag = pN2-variable, user must specify P units')
        end
        X_N2 = conditions.pN2/conditions.P_tot;
        mu_N = (G0_N2_gv(conditions.T_equilibrium, conditions.P_tot, X_N2, conditions.P_units)/2);
        mu_Ga = G0_GaN_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units) - mu_N;
        conditions.pH2 = 0;
        conditions.pNH3 = 0;
        conditions.pGa_vap = 0;

    elseif strcmp(conditions.mu_conditions_flag, 'N2-1atm')
        conditions.pN2 = 1;
        conditions.P_tot = 1;
        conditions.X_N2 = conditions.pN2/conditions.P_tot;

        mu_N = (G0_N2_gv(conditions.T_equilibrium, conditions.P_tot, conditions.X_N2, conditions.P_units)/2);
        mu_Ga = G0_GaN_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units) - mu_N;
        conditions.pH2 = 0;
        conditions.pNH3 = 0;
        conditions.pGa_vap = 0;

    elseif strcmp(conditions.mu_conditions_flag, 'air-1atm')
        conditions.pN2 = 0.8;
        conditions.P_tot = 1;
        X_N2 = conditions.pN2/conditions.P_tot;
        mu_N = (G0_N2_gv(conditions.T_equilibrium, conditions.P_tot, X_N2, conditions.P_units)/2);
        mu_Ga = G0_GaN_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units) - mu_N;
        conditions.pH2 = 0;
        conditions.pNH3 = 0;
        conditions.pGa_vap = 0;

    elseif strcmp(conditions.mu_conditions_flag, 'N2 determined from NH3/H2 Equilibrium')

        X_H2 = conditions.pH2/conditions.P_tot;
        X_NH3 = conditions.pNH3/conditions.P_tot;
        
        mu_N = 0.5*((2*G0_NH3_gv(conditions.T_equilibrium, conditions.P_tot, X_NH3, conditions.P_units)) - (3*G0_H2_gv(conditions.T_equilibrium, conditions.P_tot, X_H2, conditions.P_units)));
        conditions.pN2 = conditions.P_tot * exp((-2*mu_N)/conditions.kBT_equilbrium);
        mu_Ga = G0_GaN_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units) - mu_N;
        conditions.pGa_vap = conditions.P_tot * exp(-mu_Ga/conditions.kBT_equilbrium);

    


    else
        error('using T dependent chemical potentials, conditions.mu_conditions_flag must be one of the coded scenarios  ')
    end  % end main if elseif else structure


    % 
     
     %% these two lines make sure that all values are scalars, and that the total pressure is >= to the sum of reactive (not inert gas) species
    if ~isscalar(conditions.P_tot) || ~isscalar(conditions.pGa_vap) || ~isscalar(conditions.pN2) || ~isscalar(conditions.pH2) || ~isscalar(conditions.pNH3)     % check that each pressure is a scalar
        error('all pressures must scalar values  - make a loop over pressure in the script calling this if multiple P values are desired')
    end

    if conditions.P_tot<(conditions.pGa_vap + conditions.pN2 + conditions.pH2 + conditions.pNH3)   % check that the total pressure (whcih could contain inert gases) is greater than or equal to sum of the specified species
        error('the total pressure must be >= sum of reactive (non inert) species pressures')
    end



end

end   % end main function





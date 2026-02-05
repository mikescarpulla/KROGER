%% Function that returns the Ga and N chemical potentials in equilibrium with GaN for a vector of temperatures

function [mu_Ga, mu_N] = Ga_N2_GaN_chem_potentials_from_conditions(conditions)

%constants
kB_eV = 8.617333262e-5;

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


    conditions.pH2 = 0;
    conditions.pNH3 = 0;
    conditions.pGa_vap = 0;

elseif strcmp(conditions.T_dep_matrix_mu_flag,'On')   %% T dependent chemical potentials from thermochemistry
    disp('using T-dependent chemical potentials')

    if strcmp(conditions.mu_conditions_flag, 'GaN_touching_Ga(l,s)')
        mu_Ga = G0_Ga_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units);
        mu_N = G0_GaN_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units) - mu_Ga;
        conditions.pN2 = 0;
        conditions.pH2 = 0;
        conditions.pNH3 = 0;
        conditions.pGa_vap = 0;

    elseif strcmp(conditions.mu_conditions_flag, 'Ga_equilibrium_vapor_pressure_over_Ga_at_T')

        if (~exist('conditions.P_ref') || ~exist('conditions.P_units')) || (isempty(conditions.P_units) || isempty(conditions.P_ref))
            conditions.P_units = 'atm';  % default - give it in atm
            conditions.P_ref = 1;
        elseif strcmp(conditions.P_units,'Torr') || conditions.P_ref == 760
            conditions.P_ref = 760;
            conditions.P_units = 'Torr';
        elseif strcmp(conditions.P_units,'Bar')
            conditions.P_ref = 1;
            conditions.P_units = 'Bar';
        elseif strcmp(conditions.P_units,'Pa')|| conditions.P_ref == 1e5
            conditions.P_ref = 1e5;
            conditions.P_units = 'Pa';
        else
            error('Units of pressure must be atm, Torr, Pa, or Bar')
        end

        G0_ls = G0_Ga_ls(conditions.T_equilibrium, conditions.P_ref, 1, conditions.P_units);
        G0_vap = G0_Ga_gv(conditions.T_equilibrium, conditions.P_ref, 1, conditions.P_units);
        mu_Ga = G0_vap - G0_ls;    % mu = kBTln(vapor pressure), vapor pressure is pref*exp(-(G0 difference)/kBT) so jsut go straight to the point
        mu_N = G0_GaN_ls(conditions.T_equilibrium, conditions.P_ref, 1, conditions.P_units) - mu_Ga;

        conditions.pGa_vap = Ga_equilib_vap_pres(conditions.T_equilibrium, conditions.P_units);
        conditions.pN2 = 0;
        conditions.pH2 = 0;
        conditions.pNH3 = 0;
        conditions.P_tot = conditions.pGa_vap;

    elseif strcmp(conditions.mu_conditions_flag, 'pN2-variable') %equilibrium with pN2 set to arbitrary value. You must also set pN2.
        if isempty(conditions.pN2)
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

    elseif strcmp(conditions.mu_conditions_flag, 'NH3-H2-N2 Equilibrium')
        % rxn: N2 + 3H2 = 2NH3
        % assumption: we maintain a constant total pressure and supply
        % input gases with molar flow rates n_inert, n_NH3, n_H2, and n_N2
        % These input gases react completely and form the equilibrium mixture of H2, N2,
        % and NH3 (and inert) at total pressure Ptot.
        % inputs:
        % conditions.n_init_NH3
        % conditions.n_init_N2
        % conditions.n_init_H2
        % conditions.n_init_inert
        % conditions.T_equilibrium
        % conditions.P_tot


        conditions.NH3_N2_H2_deltaG = 2*G0_NH3_gv(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units) - G0_N2_gv(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units) - 3*G0_H2_gv(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units);  %this is at standard pressure
        % Keq = exp(-deltaG./(kB_eV*conditions.T_equilibrium));  

        N = numel(conditions.T_equilibrium);
        n_tot = zeros(N,1);
        conditions.n_final_NH3 = zeros(N,1);
        conditions.n_final_N2 = zeros(N,1);
        conditions.n_final_H2 = zeros(N,1);
        conditions.kBT_equilibrium = kB_eV*conditions.T_equilibrium;

x = -1:0.01:1;

        for i=1:N
            % conditions.Keq = Keq(i);
            delta_G = conditions.NH3_N2_H2_deltaG(i);
            kBT = conditions.kBT_equilibrium(i);
            rxn_n = fzero(@equilib_gas_mix, 0.5);
            conditions.n_final_NH3(i) = conditions.n_init_NH3 - 2*rxn_n;
            conditions.n_final_N2(i) = conditions.n_init_N2 + 2*rxn_n;
            conditions.n_final_H2(i) = conditions.n_init_H2 + 3*rxn_n;
        end

        % mole fractions in final mixture.  Ptot*X = partial pressure
        n_tot =  conditions.n_final_H2 + conditions.n_final_N2 + conditions.n_final_NH3 + conditions.n_inert;
        conditions.X_NH3 = conditions.n_final_NH3/n_tot;
        conditions.X_N2 = conditions.n_final_N2/n_tot;
        conditions.X_H2 = conditions.n_final_H2/n_tot;
        conditions.X_inert = conditions.n_inert/n_tot;
        mu_N = 0.5*G0_N2_gv(conditions.T_equilibrium, conditions.P_tot, conditions.X_N2, conditions.P_units);
        mu_Ga = G0_GaN_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units) - mu_N;

    else
        error('using T dependent chemical potentials, conditions.mu_conditions_flag must be one of the coded scenarios  ')
    end  % end main if elseif else structure


    if conditions.P_tot<(conditions.pGa_vap + conditions.pN2 + conditions.pH2 + conditions.pNH3)   % check that the total pressure (which could contain inert gases) is greater than or equal to sum of the specified species
        error('the total pressure must be >= sum of reactive (non inert) species pressures')
    end

end % T dep mu or mu const if structure


    function [Y] = equilib_gas_mix(x)
        Y = -kBT*( 2*log(conditions.n_init_NH3 - 2*x) - log(conditions.n_init_N2 + x) - 3*log(conditions.n_init_H2 + 3*x)) - delta_G;
    end




end   % end main function





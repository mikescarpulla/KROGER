%% Function that returns the Cd and Te chemical potentials in equilibrium with CdTe for a vector of temperatures

function [mu_Cd, mu_Te] = Cd_Te_CdTe_chem_potentials_from_conditions(conditions)

disp('Assumption is that the specified chemical potentials are continuously enforced at constant P (no products build up or reactants deplete). This is opposed to setting these initial conditions in a closed box and letting them react.  So user make sure specified conditions are compatible and undrstood.')

if strcmp(conditions.T_dep_matrix_mu_flag,'Off')   %T independent mu values
    disp('using T-independent chemical potentials')

    if strcmp(conditions.mu_conditions_flag,'Te-rich')
        mu_Te = 0;
        mu_Cd = conditions.mu_CdTe;

    elseif strcmp(conditions.mu_conditions_flag,'Cd-rich')   %% for now main one should be O-rich. Next the immediate next phase will be in place of "Ga"-rich
        mu_Cd = 0;
        mu_Te = conditions.mu_CdTe;
    else
        error('using T-independent mu values, conditions.mu_conditions_flag should be Cd-rich or Te-rich')
    end

    % user needs to commit these to the conditions variable in the main script
    mu_Cd = ones(size(conditions.T_equilibrium,2),1)*mu_Cd;
    mu_Te = ones(size(conditions.T_equilibrium,2),1)*mu_Te;



elseif strcmp(conditions.T_dep_matrix_mu_flag,'On')   %% T dependent chemical potentials from thermochemistry
    disp('using T-dependent chemical potentials')

    if strcmp(conditions.mu_conditions_flag,'CdTe_touching_Cd(l,s)')
        mu_Cd = G0_Cd_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units);
        mu_Te = G0_CdTe_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units) - mu_Cd;

    elseif strcmp(conditions.mu_conditions_flag,'CdTe_touching_Te(l,s)')
        mu_Te = G0_Te_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units);
        mu_Cd = G0_CdTe_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units) - mu_Te;

    elseif strcmp(conditions.mu_conditions_flag,'Cd_vapor_over_Cd_at_T_Cd')
        G_Cd_condensed = G0_Cd_ls(conditions.T_Cd, conditions.P_tot, 1, conditions.P_units);
        G_Cd_vap = G0_Cd_gv(conditions.T_Cd, conditions.P_tot, 1, conditions.P_units);
        mu_Cd = conditions.P_ref*exp(-(G_Cd_vap-G_Cd_condensed)/(kB_eV*conditions.T_Cd));  % treating Cd vapor as ideal gas at the T_Cd
        mu_Te = G0_CdTe_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units) - mu_Cd;

    elseif strcmp(conditions.mu_conditions_flag,'Te2_vapor_over_Te_at_T_Te')
        % here we approximate the vapor as being only Te2
        % dG_Te2 = 1/2*G0_Te2_gv(conditions.T_Te, conditions.P_ref, 1, conditions.P_units) - G0_Te_ls(conditions.T_Te, conditions.P_ref, 1, conditions.P_units);
        % p_Te2 = conditions.P_ref*exp(-2*(dG_Te2)/conditions.T_Te);  % vapor pressure of Te2 over Te at the Te temperature
        % mu_Te = 0.5*kB_eV*conditions.T_Te*log(p_Te2/conditions.P_ref);
        mu_Te = G0_Te_ls(T,1,1,P_units) - G0_Te2_gv(T,conditions.P_tot,1,P_units)/2;  % this avoids computing pVap then taking ln of it - the 3 lines above do it that way
        mu_Cd = G0_CdTe_ls(conditions.T_equilibrium, conditions.P_tot, 1, conditions.P_units) - mu_Te;

  
    % elseif strcmp(conditions.mu_conditions_flag,'CdTe_congruent_vapor_at_Tsource')  %as for CSS - Tequilib is the sample T while the vapor source could be higher 



    else
        error('using T dependent chemical potentials, conditions.mu_conditions_flag must be one of the coded scenarios  ')
    end  % end main if elseif else structure


    % 
    % 
    % %% these two lines make sure that all values are scalars, and that the total pressure is >= to the sum of reactive (not inert gas) species
    % if ~isscalar(conditions.P_tot) || ~isscalar(conditions.pGa_vap) || ~isscalar(conditions.pO2) || ~isscalar(conditions.pGaO_vap) || ~isscalar(conditions.pGa2O_vap)     % check that each pressure is a scalar
    %     error('all pressures must scalar values  - make a loop over pressure in the script calling this if multiple P values are desired')
    % end
    % 
    % if conditions.P_tot<(conditions.pGa_vap + conditions.pO2 + conditions.pGaO_vap + conditions.pGa2O_vap)   % check that the total pressure (whcih could contain inert gases) is greater than or equal to sum of the specified species
    %     error('the total pressure must be >= sum of reactive (non inert) species pressures')
    % end



end

end   % end main function





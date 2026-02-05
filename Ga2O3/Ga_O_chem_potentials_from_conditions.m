function [mu_Ga, mu_O] = Ga_O_chem_potentials_from_conditions(conditions)
% % function [mu_Ga, mu_O] = Ga_O_chem_potentials_from_conditions(conditions, plot11_yes_no)

% this function handles setting the Ga and O chemical potentials for the
% host material Ga2O3.  It takes in the whole conditions structure and
% returns the mu_Ga and mu_O vectors for the requested condiitons.
% enter 'y' or 'yes' to make plot 11
% Note the conditions variable has to contain some but not all of its full suite of conditions - onluy those sub-parts called herein are neededR

% handle chemical potentials of O and Ga

if strcmp(conditions.T_dep_mu_flag,'Off')   %T independent mu values
    % %%% Constant chemical potentials
    conditions.mu_Ga2O3 = -10.26;        % formation energy of Ga2O3 calc by Joel

    disp('using T-independent chemical potentials')

    if strcmp(conditions.mu_conditions_flag,'O-rich')
        mu_O = 0;
        mu_Ga = conditions.mu_Ga2O3/2;

    elseif strcmp(conditions.mu_conditions_flag,'Ga-rich')
        mu_Ga = 0;
        mu_O = conditions.mu_Ga2O3/3;

    else
        error('using T-independent mu values, conditions.mu_conditions_flag should be O-rich or Ga-rich')
    end

    % user will need to commit these to the conditions variable in the main
    % script
    mu_Ga = ones(size(conditions.T_equilibrium,2),1)*mu_Ga;
    mu_O = ones(size(conditions.T_equilibrium,2),1)*mu_O;


elseif strcmp(conditions.T_dep_mu_flag,'On')   %% T dependent chemical potentials from Zinkevich & Adlinger paper for the O and Ga

    disp('using T-dependent chemical potentials')

    if strcmp(conditions.mu_conditions_flag,'pO2-variable')   % equilibrium with pO2 set to arbitrary value.  you must also set pO2.
        conditions.pGa_vap = 0;
        conditions.pGa2O_vap = 0;
        conditions.pGaO_vap = 0;
        if isempty(conditions.pO2)
            error('if using conditions.mu_conditions_flag = pO2-variable, user must specify a pO2')
        end

    elseif strcmp(conditions.mu_conditions_flag,'pGa_vap-variable')   % equilibrium with Ga_vapor.  Must also specify a pGa_vap
        conditions.pO2 = 0;
        conditions.pGa2O_vap = 0;
        conditions.pGaO_vap = 0;
        if isempty(conditions.pGa_vap)
            error('if using conditions.mu_conditions_flag = pGa_vap-variable, user must specify a pGa_vap')
        end

    elseif strcmp(conditions.mu_conditions_flag,'O2-1atm')   % equilibrium with O2 at pO2 = 1 atm
        conditions.pGa_vap = 0;
        conditions.pGa2O_vap = 0;
        conditions.pGaO_vap = 0;
        conditions.P_units = 'atm';  % this scenario overrides any earlier setting of pressure units.
        conditions.pO2 = 1;

    elseif strcmp(conditions.mu_conditions_flag,'Ga_equilib_vap_pressure')   % equilibrium vapor pressure of Ga ovr Ga
        conditions.pGa2O_vap = 0;
        conditions.pGaO_vap = 0;
        conditions.pO2 = 0;
        conditions.pGa_vap = Ga_equilib_vap_pres(conditions.T_equilibrium, conditions.P_units);

    elseif strcmp(conditions.mu_conditions_flag,'Ga_ls-rich')   % equilibrium of Ga2O3 in direct contact with pure liquid or solid Ga
        conditions.pGa_vap = 0;
        conditions.pGa2O_vap = 0;
        conditions.pGaO_vap = 0;
        conditions.pO2 = 0;

    elseif strcmp(conditions.mu_conditions_flag,'pGa2O_vap-variable')   % equilibrium of Ga2O3 in contact with a specified p_Ga2O
        conditions.pGa_vap= 0;
        conditions.pGaO_vap = 0;
        conditions.pO2 = 0;
        if isempty(conditions.pGa2O_vap)
            error('if using conditions.mu_conditions_flag = pGa2O_vap-variable, user must specify a pGa2O')
        end



    elseif strcmp(conditions.mu_conditions_flag,'pGaO_vap-variable')   % equilibrium of Ga2O3 in contact with a specified p_GaO
        conditions.pGa_vap= 0;
        conditions.pGa2O_vap = 0;
        conditions.pO2 = 0;
        if isempty(conditions.pGaO_vap)
            error('if using conditions.mu_conditions_flag = pGa2O_vap-variable, user must specify a pGaO')
        end


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

    end

    % given the scenario and inputs, calculate the Ga and O chemical potentials
    % [mu_Ga, mu_O] = Ga2O3_chemical_potentials(conditions.T_equilibrium, conditions.P_tot, conditions.pGa_vap, conditions.pO2, conditions.pGa2O, conditions.pGaO, conditions.P_units, conditions.mu_conditions_flag);
    [mu_Ga, mu_O] = Ga2O3_chemical_potentials(conditions);


else
    error('using T dependent chemical potentials, conditions.mu_conditions_flag must be Ga-rich, O-rich, Ga2O, O-variable, or manual until other things are implemented.  ')
end  % end main if elseif else structure

% % debugging plot
% % if strcmp(plot11_yes_no,'yes') ||strcmp(plot11_yes_no,'y')
% %
% % % make a plot of chem potentials
% %     [mu_Ga_O2_1atm, mu_O_O2_1atm] = Ga2O3_chemical_potentials(conditions.T_equilibrium, 1, 0, 1, 0, 0, 'atm','O2-1atm');   % this is as rich as O2 can get at 1 atm - can increase pO2 and pTot to go higher
% %     [mu_Ga_pureGa, mu_O_pureGa] = Ga2O3_chemical_potentials(conditions.T_equilibrium, 1, 0, 0, 0, 0, 'atm','Ga_ls-rich');   % this is as rich as Ga can go when in direct contact with Ga2O3
% %     [mu_Ga_Ga_vap_press, mu_O_Ga_vap_press] = Ga2O3_chemical_potentials(conditions.T_equilibrium, 1, 0, 0, 0, 0, 'atm','Ga_equilib_vap_pressure');   % chem potentials when Ga(vap) in quilibrium with Ga(l) or Ga(s) and Ga2O3
% %
% % error('where is this called from?')
% %
% %     figure(11)
% %     clf
% %     hold on
% %     f11line(1) = plot(conditions.T_equilibrium,mu_O,'r-','DisplayName','\mu_{O}');
% %     f11line(2) = plot(conditions.T_equilibrium,mu_Ga,'r--','DisplayName','\mu_{Ga}');
% %     f11line(3) = plot(conditions.T_equilibrium,mu_O_O2_1atm,'k-','DisplayName','\mu_{O} 1 atm O2');
% %     f11line(4) = plot(conditions.T_equilibrium,mu_Ga_O2_1atm,'k--','DisplayName','\mu_{Ga} 1 atm O2');
% %     f11line(5) = plot(conditions.T_equilibrium,mu_O_pureGa,'b-','DisplayName','\mu_{O} direct contact with pure Ga');
% %     f11line(6) = plot(conditions.T_equilibrium,mu_Ga_pureGa,'b--','DisplayName','\mu_{Ga} direct contact with pure Ga');
% %     f11line(7) = plot(conditions.T_equilibrium,mu_O_Ga_vap_press,'m-','DisplayName','\mu_{O} Ga equilibrium vapor');
% %     f11line(8) = plot(conditions.T_equilibrium,mu_Ga_Ga_vap_press,'m--','DisplayName','\mu_{Ga} Ga equilibrium vapor');
% %     %         title(strcat('\mu_O and \mu_{Ga} vs T for p_{O2}=',num2str(conditions.pO2)))
% %     title('\mu_O and \mu_{Ga} vs T')
% %     xlabel('T (K)')
% %     ylabel('Chemical Potential (eV)')
% %     datacursormode.Enable='on';   %turn on the datatip to give name of defect on mouse hover
% %     for i = 1:8
% %         f11line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({f11line(i).DisplayName},size(f11line(i).XData)));
% %     end
% %     legend('\mu_{O}','\mu_{Ga}', '\mu_{O} 1 atm O_{2}','\mu_{Ga} 1 atm O2','\mu_{O} contact with pure Ga','\mu_{Ga} contact with pure Ga','\mu_{O} Ga equilib vapor','\mu_{Ga} Ga equilib vapor')
% %     clear datacursormode f11line i    %%%%%% clean up after this plot

end

function [new_concentrations] = defect_kinetics_during_cooling(current_concentrations, material, defects, conditions)



set up things here like timesteps based on the deltaT and the cooling rate.  
Loop over the T values one at a time.  Output should be a list of all concentrations at each temperature, so we can follow the trajectories of all defects vs cooling

call the ODE solver here - tell it to use @compute_dCdt as the function.  
evolve the concentrations to the next T, then do it again all the way down to the stop temperature





function [dCdt] = compute_dCdt(cs_concentrations)

% this is a calc that works for one T at a time.
% to use the built in ODE solvers, we need a function that takes in only
% the curent concentrations and outputs only the dC/dt derivatives.  

## these are here for testing, they should be passed from the conditions variable instead
kB = 8.617e-5;   % eV/K
eps0 = 8.85e-14;  % F/m
kBT = kB*T;    %eV

%compute all the current diffusion coefficients for this temperature

% loop over all chargestates and compute the diffusion prefactor.  For
% native defects and interstitials, this is jsut Do.  For things that
% diffuse mediated by some native defect, we take the standard Do and
% multiply it by the current concentration of that defect.

% compute the normal diffusion constants from Arrhenius law
D = defects.cs_D0 .* exp(-defects.cs_Ea/kBT);
mediated_diffusion_indices = find(defects.cs_diffusion_mediator_index>0);  % get the indices of the defects whose diffusion is mediated by a native defect
for i = mediated_diffusion_indices  %fix just the mediated ones in a for loop
    D(i) = D(i) * defect_concentrations(defects.cs_diffusion_mediator_index(i));  % for the defects that have mediated diffusion, multiply by the concentration of their mediators
end

dCdt = zeros(defects.num_chargestates,1);

% compute the current, actual activity coefficients of all chargestates from
% current conc and their site density prefactors  activity = #/cm3 /
% prefactor.  So this should be a dimensionless number, since Keq =
% dimensionless also.
cs_activities = cs_concentrations./defects.cs_tot_prefactorfactor;


num_rxns = size(reactions,2);
G_form_equilib = defects.cs_Ho + Ef*defects.cs_q;  % Ef should come from the equilibrium_dark_sol
G_form_actual = compute this from the current concentrations and the prefactors for each chargestate.  We can compare Gactual to G_equilib for each defect to predict if its number should be going up or down by reactions.
% then we check at the end after the time-dependent solver to see if all
% defects concentrations are headed in the right directions or if things
% look like a numerical problem that went crazy.

for rxn_index = 1:num_rxns
    stoich_coeffs = reactions(:,rxn_index);   % these are signed at this point, neg for reactants
    reactants_indices = find(stoich_coeffs<0);
    p1_index = find(stoich_coeffs>0); % this line should either find one positive number or none, so p1_index will be scalar or empty

    if numel(reactants_indices)==2
        r1_index = reactants_indices(1); %this is now the index of reactant 1 in the list of chargestates
        r2_index = reactants_indices(2);
        r1_stoich = -1*stoich_coeffs(r1_index);
        r2_stoich = -1*stoich_coeffs(r2_index);
        r1_activity = cs_activities(r1_index);
        r2_activity = cs_activities(r2_index);
    elseif numel(reactants_indices)==1  %is the reactions only lists one reactant, it must be reacting with itself to form a pair
        if stoich_coeffs(reactants_indices)==2  %if the thing reacts with itself (like a pair of vacancies)
            r1_index = reactants_indices(1);
            r2_index = reactants_indices(1);
            r1_stoich = -1*stoich_coeffs(r1_index);
            r2_stoich = -1*stoich_coeffs(r2_index);
            r1_activity = cs_activities(r1_index);
            r2_activity = cs_activities(r2_index);
        end
    else
        error('Reactions must be written so there is one or more chargestate as a reactant, and as binary reactions so the numbers in the reaction list have to be 0, -1, or -2')
    end

    % now we have identified our reactants.  Usually two oppositely charged
    % things will attract eachother.  But we can also have like charges
    % come together, or neutrals come together
    radius(rxn_index) = abs(defects.cs_q(reactants_indices(1))*defects.cs_q(reactants_indices(2))) /(4*pi()*eps0*material.eps_rel*kBT);    radius = -(q1 q2) / (4 pi eps eps0 kBT)  % onsager escape radius for thermal energy compared to Coulomb attraction (note if chargestates have like charge, we assume rxn doesnt happen so the reaction would not be listed in the first place)
    if radius(rxn_index)==0  %if the two things are neutral
        radius(rxn_index)=min_radius;  ### define this somewhere in the material variable
    elseif radius(rxn_index)<0   %somehow if eps or some other number is accidentally negative
        radius(rxn_index)=min_radius;
    end


    if numel(p1_index)==1  %normal case
        p1_stoich = stoich_coeffs(p1_index);
        p1_activity = cs_activities(p1_index);
    elseif isempty(product_index)   % case where reactants annihilate eachother like interstitial plus vacacny
        p1_stoich = 1;
        p1_activity = 1;
    else
        error('Reactions must be written so there is one product or no products for annihilation reactions')
    end

    % compute the deltaG for the reaction and the Keq
    delta_G_rxn = G_form_equilib(r1_index) - G_form_equilib(r1_index) - G_form_equilib(r2_index);
    Keq = exp(-delta_G_rxn/kBT);

    kfwd = radius * ( D(r1_index) + D(r2_index)) ;   % cm * cm2/s = cm3/s
    kback = kfwd/Keq;  % impose detailed balance

    Q_num = p1_activity^p1_stoich;
    Q_denom = (r1_activity^r1_stoich) * (r2_activity^r2_stoich);
    Q_rxn = Q_num/Q_denom;  % this is the current value of the reaction quotient.  The reaction proceeds to right if it is >Keq and to left if it is <Keq
    rxn_direction = sign(log(Q_rxn/Keq));  % if this is negative, it means reaction goes backwards.  If positive, reaction should go forwards.

    rf = kfwd * Q_denom;
    rb = kback * Q_num;

    net_rxn_rate(rxn_index) = (rf-rb)/p1_stoich;  % this is the net rate of this reaction, normalized to the stoichiometric coefficient of the product

    % compute the rates reactants disappear and products appear.  each time
    % through the loop over reactions we may add or subtract from the net
    % dC/dt for each chargestate
    dCdt(r1_index) = dCdt(r1_index) + net_rxn_rate(rxn_index) * r1_stoich;
    dCdt(r2_index) = dCdt(r2_index) + net_rxn_rate(rxn_index) * r2_stoich;
    dCdt(p1_index) = dCdt(p1_index) + net_rxn_rate(rxn_index) * p1_stoich;

end  %end looping over the reactions

end %end the dC/dt = F( C ) function for the ODE solver




end %end the main function
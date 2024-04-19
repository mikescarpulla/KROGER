% version of defect equilibrium that does calc for just one temperature,
% not a vector of them (since each needs different treatment in GQ).
% Anticipated use is to send in gq_defects instead of main defects
% structure


function [gq_equilib_dark_sol] = defect_gq_equilibrium_dark(gq_Tnow, material, gq_defects, conditions)

kB = 8.6173324e-5; 
N_charge_states_out = zeros(1,gq_defects.numchargestates);   % in the final output, each charge state has a column, each row is a different T. Need this different name for the holder matrix because of scope of n_defects in some of the functions
charge_bal_err_out = zeros(1,1);
kBT = gq_Tnow*kB;    % this is just like the regular equilibrium dark calc but this temperature overrides the vector of T's (conditions.T) usually used.  Thus this function will only pop out one line of resutls for Tnow.   

[guess] = EF_guess(kBT,material.Eg);   % get a guess for EF close to minimum using trial and error
%%% find root of charge bal equation
[EF_out, ~,~,~] = fzero(@charge_bal,guess);
%     [EF_out, fval, exitflag, output] = fzero(@charge_bal,guess);
[n_out, p_out, N_charge_states_out] = concentrations(EF_out,kBT);  % compute the defect and carrier concentrations from the Ef
[charge_bal_err_out] = charge_bal(EF_out);   % compute the net charge just to check final answers
%%%%%%%%%%%%%% end of main calculations  %%%%%%%%%%%%%%%%%%


% build up the solution structure
N_defects_out = zeros(1,gq_defects.numdefects);
for i=1:gq_defects.numdefects
    index = gq_defects.ID==i;
    N_defects_out(1,i) = sum(N_charge_states_out(1,index),2);
end


gq_equilib_dark_sol.Nd = material.Nd;
gq_equilib_dark_sol.Na = material.Na;
gq_equilib_dark_sol.Efn = EF_out;   % planning ahead for light calcs
gq_equilib_dark_sol.Efp = EF_out;
gq_equilib_dark_sol.n = n_out;
gq_equilib_dark_sol.p = p_out;
gq_equilib_dark_sol.chargestates = N_charge_states_out;
gq_equilib_dark_sol.defects = N_defects_out;
gq_equilib_dark_sol.charge_bal_err = charge_bal_err_out;
gq_equilib_dark_sol.defectnames = gq_defects.defectnames;
gq_equilib_dark_sol.chargestatenames = gq_defects.chargestatenames;


%%%%%%%%% end of main function, end kept open to make subroutines nested  %%%%%%%%%%%%%%%  



%%%%%%%%%%%%%%% nested subroutines %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% function that uses grid search to get Ef close to the charge
%%%% balance solution.  It will output two values for Ef that bracket the solution unless something is strange about the charge_bal vs Ef (i.e. its not monotonic   %%%%%%%%%%%%
% could change this to just guess based on net doping rather than grid
% seach
function [guess] = EF_guess(kBT,Egap)   % keep kBT and Eg as local variables in this function

    EF_int = kBT/2;   % this guarantees you can't miss the solution whcih should be thus within kB/2 of the guess
%     EF_grid = (-5*kBT):EF_int:(ceil(Egap/EF_int)*EF_int + 5*kBT);  % this makes a grid to check over the range -5kBT to Eg+5kBT 
    EF_grid = 2.5*kBT:EF_int:(ceil((Egap-2.5*kBT)/EF_int))*EF_int;  % this makes a grid to check over the range for which Boltzmann approx is OK
%     EF_grid = Egap/2:EF_int:(ceil(Egap/EF_int)*EF_int);  % this makes a grid to check over the range Eg/2 to Eg
    nn = size(EF_grid,2);
    errs = zeros(1,nn);

    for i=1:nn
        guess = EF_grid(i);
        errs(i)=charge_bal(guess);
        % errs(i) = square_err(guess);
    end
    
    edge_index = find(diff(sign(errs))~=0);   %  this finds the rising/falling edge where error changes sign
    
    if size(edge_index)==[1 1]
        min_index = [edge_index edge_index+1];   % so we are finding the two guesses that bracket the solution one + and one - 
    elseif size(edge_index)~=[1 1]
        disp('waring: charge balance error may not be monotonic')
%     else
%         error('something strange about charge balance error vs Ef - solutions may not be valid')
    end
    
    
    
%     min_index = find(errs==min(errs),1,'first');  %find the grid value which is closest to real answer and pass this to the main calculation  -use this for case of function minimization
%     min_index  = find(abs(errs)==min(abs(errs)),1,'first');  % find the grid value that gives the closest value to 0 for charge bal - use this version for root/zero finding method.  
%     [~,min_index] = min(abs(errs),[],1,'linear');  % find the grid value that gives the closest value to 0 for charge bal - use this version for root/zero finding method.  
%     disp(min_index)
    guess = EF_grid(min_index);
%      
% %%%% plot the charge balance error function vs EF position if you like
% figure(1)
% clf
% plot(EF_grid,errs)
% hold on
% % plot(EF_grid(min_index),errs(min_index),'ro')  %plot the EF_guess 


end
%%%% end EF_guess 





%%%% function that computes the net charge.   N_chargestates is a local variable during solution process 
function charge_bal = charge_bal(EF)
    [n, p, N_chargestates] = concentrations(EF,kBT);
    charge_bal = (sum( gq_defects.charge.* N_chargestates') + p - n + material.Nd - material.Na);    

end
%%%% end charge_bal   %%%%%%%%%%%%%








 
%%% version without frozen defects
%%%%  function that calculates the numbers of each charge state of each defect plut electrons and holes  %%%%%%%%%%%%%%%%%%%
%% note here in future the 1/(1+exp(-dH1)+exp(-dH2)+exp(-dH3)+...) site blocking could be implemented.  Would need to use defects.sites.  Then loop over each site in the lattice and find all of the charge states that try to occupy it, then construct the FD-like function for that site.  Then use that.  Might be faster to compute exp(-dH/kBT) for each charge state first, then sum 1+the rights ones for each distribution function.  
function [n, p, N_chargestates] = concentrations(EF,kBT)
dH = gq_defects.dHo - gq_defects.dm*conditions.mu + gq_defects.charge*EF;

N_chargestates = (gq_defects.pre .* exp(-dH/kBT))';   % Boltzmann dilute defects approximation here for all charge states of all defects.  
[n,p] = carrier_conc(EF, kBT);

end
%%%% end concentrations



% % %%%%  function that calculates the numbers of each type of defects, electrons, and holes  %%%%%%%%%%%%%%%%%%%
% % %% this one incudes frozen defects
% % %% note here in future the 1/(1+exp(-dH1)+exp(-dH2)+exp(-dH3)+...) site blocking could be implemented.  Would need to use defects.sites.  Then loop over each site in the lattice and find all of the charge states that try to occupy it, then construct the FD-like function for that site.  Then use that.  Might be faster to compute exp(-dH/kBT) for each charge state first, then sum 1+the rights ones for each distribution function.  
% % function [n, p, N_chargestates] = concentrations(EF,kBT)
% % 
% %     % calc all the chargestates according to usual formation enthalpy calc
% %     dH = gq_defects.dHo - gq_defects.dm*conditions.mu + gq_defects.charge*EF;
% %     N_chargestates = (gq_defects.pre.*exp(-dH/kBT))';   % Boltzmann dilute defects approximation here for all charge states of all defects.
% %     [n,p] = carrier_conc(EF, kBT);
% % 
% %     if sum(gq_defects.gq_isfrozen_allT_defects)~=0
% %         % now handle the defects that ARE frozen - just overwrite the values calculated under usual open system in lines above.  Have to do a partition function calc for them.  The sum of all chargestates within a defect will be set by the specified concentration for that defect
% %         % calculate Boltzmann factors for the frozen chargestates using formation enthalpy without the chem potential term
% %         dH_rel = zeros(gq_defects.numchargestates,1);
% %         Boltz_facs = dH_rel;
% %         dH_rel(gq_defects.gq_isfrozen_allT_chargestates_index) = gq_defects.dHo(gq_defects.gq_isfrozen_allT_chargestates_index) + gq_defects.charge(gq_defects.gq_isfrozen_allT_chargestates_index)*EF;
% %         Boltz_facs(gq_defects.gq_isfrozen_allT_chargestates_index) = exp(-dH(gq_defects.gq_isfrozen_allT_chargestates_index)/kBT); % so the unfrozen ones should be zeros, frozen ones should be finite
% % 
% %         % loop over the frozen defects and calc the concentrations
% %         Z = zeros(1,gq_defects.gq_isfrozen_allT_numfrozendefects);
% %         for i = 1:gq_defects.gq_isfrozen_allT_numfrozendefects   % loop over the frozen defects (defect.ID) - not over charge states
% %             frozen_cs_indices = gq_defects.cs_indices_lo(i):gq_defects.cs_indices_hi(i);  % vector of indices of the charge states of the ith frozen defect
% %             Z(i) = sum(Boltz_facs(frozen_cs_indices)); % Z value for the ith frozen defect (computed from Boltz factors for each charge state in that defect
% %             N_chargestates(frozen_cs_indices) = gq_defects.gq_isfrozen_allT_concdefects(i)*(Boltz_facs(frozen_cs_indices)/Z(i)) ;     % denom is a scalar.  Compute the conc of each charge state in defect i is Boltz factor/Z(i)*total conc
% %         end
% %     end
% % end
% % %%%% end concentrations










%%%% function to compute n and p from Ef  %%%%%%%%%%%%
function [n,p] = carrier_conc(EF, kBT)   % keep these variables local 
    
    etaCB = (EF - material.Ec)/kBT;
    etaVB = (EF - material.Ev)/kBT;   % this looks wrong (not symmetric compared to CB case) but it is right. Direction of integration and sign on Ef-Ev are both swapped.   
   
    
    % use just Boltzmann approx 
    n = material.Nc*exp(etaCB);   % these are right (sign swap).  Boltzmann factors should end up <1 when Ef is in gap
    p = material.Nv*exp(-etaVB);

    
    
    
%     % use just FD approximaiton
%     n = material.Nc*FD_int_approx(etaCB,1/2);
%     p = material.Nv*FD_int_approx(etaCB,1/2);
     
    
%%%% use mix of Boltzmann and FD    
%     if material.Ec-EF>3*kBT && EF>3*kBT
%         n = material.Nc*exp(etaCB);   % these are right (sign swap).  Boltzmann factors should end up <1 when Ef is in gap
%         p = material.Nv*exp(-etaVB);
%        
%     elseif material.Ec-EF<=3*kBT || EF<=3*kBT
% %         FD_order_half = Lundstrom_FD_integral_table(0.5,
% %         strcat(pname,fname);   % this is a table-lookup version of the FD
% %         function, more complicated but more accurate
% 
%         n = material.Nc*FD_int_approx(etaCB,1/2);
%         p = material.Nv*FD_int_approx(etaCB,1/2);
%         
%     else
%         error('Ef is out of bounds for the Fermi-Dirac or Boltzmann function')
%     end
        
end

%%%% end Fermi Dirac   %%%%%%%%%%%%%%%







% %%%%  function that calculates the log of the normalized square error in charge balance for a given EF, defect charges, and defect numbers  %%%%%%%%%%%%%%%%%%%
% function [charge_bal_err] = square_err(EF)
% 
% [n, p, n_defects] = concentrations(kBT, EF);
% charge_bal_err = log(( (sum( q.* n_defects') + p - n) / (sum( abs(q).* n_defects')+p +n))^2);  %using the log of the error makes the function smoother, so easier to take derivatives for the minimization function
% 
% end
% %%%% end square_err  %%%%%%%%%%%%%%%




end   %%%%%%% end main function defect equilib.  Putting the end here makes all the other functions nested in the main one and thus able to see the variables material, defects, and conditions  %%%%%%%%%%%%%%




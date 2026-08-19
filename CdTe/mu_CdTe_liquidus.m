function [Dmu_out, T_out] = mu_CdTe_liquidus(T_in, rich_flag)
% function that gives the mu_Cd and mu_Te along the Cd-Te T-x liquidus line
% Broken into Cd-rich and Te-rich sides
% T = column vector of T values
% rich_flag = "Cd-rich" or "Te-rich" depending on which side of T-x diagram
% mu_Cd, mu_Te = col vectors of mu values

%% uncomment to test the function as a script
% T_in = (1000:20:1400)';
% % rich_flag = "Cd-rich";
% rich_flag = "Te-rich";

kB = 8.6173324e-5;
Tm_Cd = 594;
Tm_Te = 722.6;
Tm_CdTe = 1365;

T_in = T_in(:);

%% compute vapor pressures of elements from thermochem data, at the T points input
dG_Cd = G_Cd_gv(T_in', 1, 1, 'atm') - G_Cd_ls(T_in', 1, 1, 'atm');
dG_Te2 = G_Te2_gv(T_in', 1, 1, 'atm') - 2*G_Te_ls(T_in', 1, 1, 'atm'); % have to account for the 2 atoms of Te vs 1 molecule of Te2
pvap_Cd = exp(-dG_Cd./(kB*T_in));
pvap_Te2 = exp(-dG_Te2./(kB*T_in));
% 
% figure(1)
% clf
% hold on
% plot(1./T_in,log10(pvap_Cd))
% plot(1./T_in,log10(pvap_Te2))
% 
% figure(2)
% clf
% hold on
% plot(T_in,log10(pvap_Cd))
% plot(T_in,log10(pvap_Te2))


%% Cd on Cd-rich side
T_inv_Cd_Cd_rich = 1/1000*[0.73241;     0.732600733;     0.733184315;     0.733917075;     0.734796926;     0.736268497;     0.742195197;     0.752237666;     0.759924207;     0.77166;     0.7829;     0.79966;     0.82017;     0.86614;     0.87992;     0.89447;     0.90869;     0.91916;     0.93085;     0.94463;     0.96051;     0.97396;     0.99148;     0.99942;     1.0993;     1.18077;     1.219512195;     1.25;     1.282051282;     1.315789474;     1.351351351;     1.388888889];
log10_p_Cd_Cd_rich = [0.09057;     0.100370545;     0.235192448;     0.33157647;     0.407148462;     0.486809953;     0.629136055;     0.742901585;     0.80062406;     0.81588;     0.80962;     0.76105;     0.68429;     0.46497;     0.39918;     0.33338;     0.25662;     0.20336;     0.1501;     0.07334;     -0.00499;     -0.07549;     -0.16478;     -0.20081;     -0.73501;     -1.16268;     -1.339836835;     -1.480855665;     -1.625444969;     -1.773790256;     -1.92609188;     -2.082566666];
p_Cd_Cd_rich = 10.^log10_p_Cd_Cd_rich;
T_Cd_Cd_rich = 1./T_inv_Cd_Cd_rich;
[T_Cd_Cd_rich, idx] = sort(T_Cd_Cd_rich(:));
p_Cd_Cd_rich = p_Cd_Cd_rich(idx);
T_inv_Cd_Cd_rich = T_inv_Cd_Cd_rich(idx);
log10_p_Cd_Cd_rich = log10_p_Cd_Cd_rich(idx);
[T_Cd_Cd_rich, idx] = unique(T_Cd_Cd_rich, 'stable');
p_Cd_Cd_rich = p_Cd_Cd_rich(idx);
T_inv_Cd_Cd_rich = T_inv_Cd_Cd_rich(idx);
log10_p_Cd_Cd_rich = log10_p_Cd_Cd_rich(idx);
clear idx
minT_Cd_Cd_rich = min(T_Cd_Cd_rich);
maxT_Cd_Cd_rich = max(T_Cd_Cd_rich);
% figure(1)
% plot(T_inv_Cd_Cd_rich, log10_p_Cd_Cd_rich,'bx');
% figure(2)
% plot(T_Cd_Cd_rich, log10_p_Cd_Cd_rich,'bx');
% 

%% Te2 on Cd-rich side
T_inv_Te2_Cd_rich = 1/1000*[0.73284;     0.73427;     0.73608;     0.73971;     0.74658;     0.75412;     0.7628;     0.7712;     0.78169;     0.78856;     0.80736];
log10_p_Te2_Cd_rich = [-2.57643;     -2.90948;     -3.21506;     -3.6185;     -4.00992;     -4.37387;     -4.67945;     -4.96099;     -5.24597;     -5.43482;     -5.87259];
p_Te2_Cd_rich = 10.^log10_p_Te2_Cd_rich;
T_Te2_Cd_rich = 1./T_inv_Te2_Cd_rich;
[T_Te2_Cd_rich, idx] = sort(T_Te2_Cd_rich(:));
p_Te2_Cd_rich = p_Te2_Cd_rich(idx);
T_inv_Te2_Cd_rich = T_inv_Te2_Cd_rich(idx);
log10_p_Te2_Cd_rich = log10_p_Te2_Cd_rich(idx);
[T_Te2_Cd_rich, idx] = unique(T_Te2_Cd_rich, 'stable');
p_Te2_Cd_rich = p_Te2_Cd_rich(idx);
T_inv_Te2_Cd_rich = T_inv_Te2_Cd_rich(idx);
log10_p_Te2_Cd_rich = log10_p_Te2_Cd_rich(idx);
clear idx
minT_Te2_Cd_rich = min(T_Te2_Cd_rich);
maxT_Te2_Cd_rich = max(T_Te2_Cd_rich);
% figure(1)
% plot(T_inv_Te2_Cd_rich, log10_p_Te2_Cd_rich,'ro');
% figure(2)
% plot(T_Te2_Cd_rich, log10_p_Te2_Cd_rich,'ro');


%% Te2 on Te-rich side
T_inv_Te2_Te_rich = 1/1000*[0.73284;     0.73389;     0.73608;     0.73751;     0.73999;     0.74257;     0.74868;     0.755;     0.76;     0.76213;     0.764;     0.77;     0.78;     0.79;     0.8;     0.80994;     0.83;     0.85;     0.88083;     0.90078;     0.95;     1;     1.05;     1.08026;     1.20011;     1.219512195;     1.25;     1.282051282;     1.315789474;     1.351351351;     1.388888889];
log10_p_Te2_Te_rich = [-2.57643;     -2.23823;     -1.95153;     -1.85368;     -1.69746;     -1.59445;     -1.38501;     -1.22184875;     -1.124938737;     -1.07943;     -1.055517328;     -0.991399828;     -0.903089987;     -0.829738285;     -0.795880017;     -0.76011;     -0.725842151;     -0.725842151;     -0.81848;     -0.89745;     -1.086186148;     -1.337242168;     -1.585026652;     -1.78329;     -2.49231;     -2.601810335;     -2.765885476;     -2.934114856;     -3.106714315;     -3.283916965;     -3.465975081];
T_Te2_Te_rich = 1./T_inv_Te2_Te_rich;
p_Te2_Te_rich = 10.^log10_p_Te2_Te_rich;
[T_Te2_Te_rich, idx] = sort(T_Te2_Te_rich(:));
p_Te2_Te_rich = p_Te2_Te_rich(idx);
T_inv_Te2_Te_rich = T_inv_Te2_Te_rich(idx);
log10_p_Te2_Te_rich = log10_p_Te2_Te_rich(idx);
[T_Te2_Te_rich, idx] = unique(T_Te2_Te_rich, 'stable');
p_Te2_Te_rich = p_Te2_Te_rich(idx);
T_inv_Te2_Te_rich = T_inv_Te2_Te_rich(idx);
log10_p_Te2_Te_rich = log10_p_Te2_Te_rich(idx);
clear idx
minT_Te2_Te_rich = min(T_Te2_Te_rich);
maxT_Te2_Te_rich = max(T_Te2_Te_rich);
% figure(1)
% plot(T_inv_Te2_Te_rich, log10_p_Te2_Te_rich,'rx');
% figure(2)
% plot(T_Te2_Te_rich, log10_p_Te2_Te_rich,'rx');

%% Cd on Te-rich side
T_inv_Cd_Te_rich = 1/1000*[1.01783;     0.99942;     0.97848;     0.96095;     0.93713;     0.9217;     0.89954;     0.87948;     0.85445;     0.83395;     0.81432;     0.79503;     0.77838;     0.75865;     0.74785;     0.74035;     0.73903;     0.73825;     0.7366;     0.73274;     0.73241];
log10_p_Cd_Te_rich = [-4.71406;     -4.501;     -4.23312;     -4.01381;     -3.70519;     -3.52034;     -3.22896;     -2.95012;     -2.60547;     -2.31409;     -2.01802;     -1.71567;     -1.41802;     -1.00759;     -0.68174;     -0.41386;     -0.34337;     -0.29481;     -0.21178;     0.00127;     0.09057];
p_Cd_Te_rich = 10.^log10_p_Cd_Te_rich;
T_Cd_Te_rich = 1./T_inv_Cd_Te_rich;
[T_Cd_Te_rich, idx] = sort(T_Cd_Te_rich(:));
p_Cd_Te_rich = p_Cd_Te_rich(idx);
T_inv_Cd_Te_rich = T_inv_Cd_Te_rich(idx);
log10_p_Cd_Te_rich = log10_p_Cd_Te_rich(idx);
[T_Cd_Te_rich, idx] = unique(T_Cd_Te_rich, 'stable');
p_Cd_Te_rich = p_Cd_Te_rich(idx);
T_inv_Cd_Te_rich = T_inv_Cd_Te_rich(idx);
log10_p_Cd_Te_rich = log10_p_Cd_Te_rich(idx);
clear idx
minT_Cd_Te_rich = min(T_Cd_Te_rich);
maxT_Cd_Te_rich = max(T_Cd_Te_rich);
% figure(1)
% plot(T_inv_Cd_Te_rich, log10_p_Cd_Te_rich,'bo');
% figure(2)
% plot(T_Cd_Te_rich, log10_p_Cd_Te_rich,'bo');

      % Note on Ptot: Greenberg's experiments were fixed volume so Ptot=pCd+pTe2.
        %  So Ptot was not 1 atm = Pref.  But ln(pCd/Ptot) + ln(ptot/pref) =
        %  ln(pCd/pref).  So the pTot cancels out.




%% fitting to interpolate the actual pCd and pTe2 vs T functions
%% These CAN NOT be used to extrapolate at all - only interpolate!!
%% use the min and maxT values for each of 4 data sets to make sure never to go outside of the interpolation range!!
p_Cd_Cd_rich_fit = @(Tq) interp1(T_Cd_Cd_rich, p_Cd_Cd_rich, Tq, 'pchip');
p_Cd_Te_rich_fit = @(Tq) interp1(T_Cd_Te_rich, p_Cd_Te_rich, Tq, 'pchip');
p_Te2_Cd_rich_fit = @(Tq) interp1(T_Te2_Cd_rich, p_Te2_Cd_rich, Tq, 'pchip');
p_Te2_Te_rich_fit = @(Tq) interp1(T_Te2_Te_rich, p_Te2_Te_rich, Tq, 'pchip');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure(1)
% plot(1./T_in, log10(p_Cd_Cd_rich_fit(T_in)),'g-')
% % plot(1./T_in, log10(p_Te2_Cd_rich_fit(T_in)),'g-')
% plot(1./T_in, log10(p_Te2_Te_rich_fit(T_in)),'g-')
% % plot(1./T_in, log10(p_Cd_Te_rich_fit(T_in)),'g-')
% xlabel("1/T (K^{-1})")
% ylabel("log10(p/atm)")
% 
% figure(2)
% plot(T_in, log10(p_Cd_Cd_rich_fit(T_in)),'g-')
% % plot(T_in, log10(p_Te2_Cd_rich_fit(T_in)),'g-')
% plot(T_in, log10(p_Te2_Te_rich_fit(T_in)),'g-')
% % plot(T_in, log10(p_Cd_Te_rich_fit(T_in)),'g-')
% xlabel("1/T (K^{-1})")
% ylabel("log10(p/atm)")

% %% Note how crazy the interpolants go if T is outside of the range for which the data exists.  We can not extrapolate at all.


%% generate the actual outputs
T_out = T_in(T_in<=Tm_CdTe);  % throw out any points above Tm_CdTe
% T_out = sort(unique(T_out));
% T_out = T_out(:);


if strcmp(rich_flag,"Cd-rich")
    % take only T values for which both the Cd and Te interpolations are valid.
    T_out = T_out(T_out >= minT_Cd_Cd_rich);
    T_out = T_out(T_out <= maxT_Cd_Cd_rich);
    % T_out = T_out(T_out>=minT_Te2_Cd_rich);
    % T_out = T_out(T_out<=maxT_Te2_Cd_rich);   

    % if the data exists for the requested temperatures, do these steps
    % % generate the pvap for Cd and Te2 at the specified T_out
    dG_Cd = G_Cd_gv(T_out,1,1,'atm') - G_Cd_ls(T_out,1,1,'atm');
    % dG_Te2 = G_Te2_gv(T_out,1,1,'atm') - 2*G_Te_ls(T_out,1,1,'atm');  %2 atoms of Te per Te2 dimer
    pvap_Cd_out = exp(-dG_Cd./(kB*T_out));
    % pvap_Te2_out = exp(-dG_Te2./(kB*T_out));
    ln_term_Cd_out = kB*T_out .* log(p_Cd_Cd_rich_fit(T_out) ./ pvap_Cd_out);
    % ln_term_Te_out = 1/2*kB*T_out .* log(p_Te2_Cd_rich_fit(T_out) ./ pvap_Te2_out);
    % 
    % %% convert to the desired delta_mu space for defect calcs - this is done by adding G0_Cd + G0_Te to both sides
    Dmu_Cd = ln_term_Cd_out + G_Cd_ls(T_out,1,1,'atm');
    % Dmu_Te = ln_term_Te_out + G_Te_ls(T_out,1,1,'atm');
    

    % in this case, the data for pTe2 doesnt go low enough in T, but good news
    % is that the pCd matches that for elemental Cd really well.  So we can
    % just use the Dmu_Cd + Dmu_Te = G_CdTe
Dmu_Te = G_CdTe_ls(T_out, 1, 1, 'atm') - Dmu_Cd;





elseif strcmp(rich_flag,"Te-rich")
    % take only T values for which both interpolations are valid
    T_out = T_out(T_out>=minT_Cd_Te_rich);
    T_out = T_out(T_out<=maxT_Cd_Te_rich);
    T_out = T_out(T_out>=minT_Te2_Te_rich);
    T_out = T_out(T_out<=maxT_Te2_Te_rich);

    % generate the pvap for Cd and Te2 at the specified T_out
    dG_Cd = G_Cd_gv(T_out,1,1,'atm') - G_Cd_ls(T_out,1,1,'atm');
    dG_Te2 = G_Te2_gv(T_out,1,1,'atm') - 2*G_Te_ls(T_out,1,1,'atm');  %2 atoms of Te per Te2 dimer
    pvap_Cd_out = exp(-dG_Cd./(kB*T_out));
    pvap_Te2_out = exp(-dG_Te2./(kB*T_out));
    ln_term_Cd_out = kB*T_out .* log(p_Cd_Te_rich_fit(T_out) ./ pvap_Cd_out);
    ln_term_Te_out = 1/2*kB*T_out .* log(p_Te2_Te_rich_fit(T_out) ./ pvap_Te2_out);

    %% convert to the desired delta_mu space for defect calcs - this is done by adding G0_Cd + G0_Te to both sides
    Dmu_Cd = ln_term_Cd_out + G_Cd_ls(T_out,1,1,'atm');
    Dmu_Te = ln_term_Te_out + G_Te_ls(T_out,1,1,'atm');
    
    % in this case, the data for pTe2 doesnt go low enough in T, but good news
    % is that the pCd matches that for elemental Cd really well.  So we can
    % just use the Dmu_Cd + Dmu_Te = G_CdTe

T_out_low = T_in(T_in < min(T_out));
Dmu_Te_low = G_Te_ls(T_out_low, 1, 1, 'atm');
Dmu_Cd_low = G_CdTe_ls(T_out_low, 1, 1, 'atm') - Dmu_Te_low;

T_out = [T_out; T_out_low];
Dmu_Te = [Dmu_Te; Dmu_Te_low];
Dmu_Cd = [Dmu_Cd; Dmu_Cd_low];



else
    error("rich_flag must be Cd-rich or Te-rich")
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% calculate the standard formation reaction delta G for CdTe from elements (in liquid states, contextually) as function of T
DG_form_CdTe_thermo = G_CdTe_ls(T_out,1,1,'atm') - G_Cd_ls(T_out,1,1,'atm')-G_Te_ls(T_out,1,1,'atm') ;



T_out = T_out';


% put the two Dmu's into 2 columns
Dmu_out = [Dmu_Cd Dmu_Te];  

% % check how the mu vlaues compare to the G_CdTe 
% DG_CdTe_out = ln_term_Cd_out + ln_term_Te_out;



% 
% figure(3)
% clf
% hold on
% plot(T_out,ln_term_Cd_out,'b-')
% plot(T_out,ln_term_Te_out,'r-')
% plot(T_out,DG_CdTe_out,'g-');
% plot(T_out,DG_form_CdTe_thermo,'k-')
% xlim([700 1400])
% ylim([-3 0])
% legend("kBT ln(p_{Cd}(T)/p_{Cd,vap}(T))", "kBT/2 ln(p_{Te2}(T)/p_{Te2,vap}(T))", "sum of ln() terms", "\DeltaG_{formation, CdTe}(T)")
% xlabel("T (K)")
% ylabel("Energy (eV)")



% 
% figure(4)
% clf
% hold on
% plot(T_out,Dmu_Cd,'b-')
% plot(T_out,Dmu_Te,'r-')
% plot(T_out,Dmu_Cd+Dmu_Te,'g-')
% plot(T_out,G_CdTe_ls(T_out,1,1,'atm'),'k')
% xlim([700 1400])
% ylim([-3 0])
% xlabel("T (K)")
% ylabel("\Delta\mu (eV)")
% legend("\Delta\mu_{Cd}","\Delta\mu_{Te}","\Delta\mu_{Cd} + \Delta\mu_{Te}", "G^o_{CdTe, thermochem}")
% 



end
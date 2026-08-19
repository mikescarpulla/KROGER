P_units = 'atm';
kB_eV = 8.617e-5;
P_ref = 1;

% Te species over elemental Te
figure(1)
clf
hold on
ylim([-20 5])
title('Vapor Pressure of Te_i (i=1..7) over Te')
xlabel('T (K)')
ylabel('log_{10}(p in atm)')



T = 300:10:1600;
% T = 722.66;  
T=T';
G_Te_condensed = G_Te_ls(T,1,1,P_units);  % this is the partial molar gibbs energy per Te in solid or liquid Te

% here we set the partial molar gibbs energies per Te equal in the two
% phases (vapor of various cluser size and the condensed phase)
% this is the change in partial molar gibbs energy per Te atom in
% solid/liquid upon vaporization.  Equations written for one atom of solid/liquid turning into 1/m molecuels of the vapor species  
dGTe1 = G_Te_gv(T,1,1,P_units) - G_Te_condensed;    
dGTe2 = 1/2*G_Te2_gv(T,1,1,P_units) - G_Te_condensed;
dGTe3 = 1/3*G_Te3_gv(T,1,1,P_units) - G_Te_condensed;
dGTe4 = 1/4*G_Te4_gv(T,1,1,P_units) - G_Te_condensed;
dGTe5 = 1/5*G_Te5_gv(T,1,1,P_units) - G_Te_condensed;
dGTe6 = 1/6*G_Te6_gv(T,1,1,P_units) - G_Te_condensed;
dGTe7 = 1/7*G_Te7_gv(T,1,1,P_units) - G_Te_condensed;

% Te above Te
% For each species we have the equation Te(s,l) = 1/m Te_m(vap) with
% Keq=exp(-dG/kBT) = (pTe_m)^(1/m).  So if you want pTe_m the pressure of
% each species, we need to multiply the dG above by m.  
p_Te1 = P_ref*exp(-dGTe1./(kB_eV*T));   
p_Te2 = P_ref*exp(-2*dGTe2./(kB_eV*T));
p_Te3 = P_ref*exp(-3*dGTe3./(kB_eV*T));
p_Te4 = P_ref*exp(-4*dGTe4./(kB_eV*T));
p_Te5 = P_ref*exp(-5*dGTe5./(kB_eV*T));
p_Te6 = P_ref*exp(-6*dGTe6./(kB_eV*T));
p_Te7 = P_ref*exp(-7*dGTe7./(kB_eV*T));
p_Tetot = p_Te1 + p_Te2 + p_Te3 + p_Te4 + p_Te5 + p_Te6 + p_Te7;

line(1) = plot(T,log10(p_Te1),'DisplayName','FactSage Te1');
line(2) = plot(T,log10(p_Te2),'ro','DisplayName','FactSage Te2');
line(3) = plot(T,log10(p_Te3),'DisplayName','FactSage Te3');
line(4) = plot(T,log10(p_Te4),'DisplayName','FactSage Te4');
line(5) = plot(T,log10(p_Te5),'DisplayName','FactSage Te5');
line(6) = plot(T,log10(p_Te6),'DisplayName','FactSage Te6');
line(7) = plot(T,log10(p_Te7),'DisplayName','FactSage Te7');
line(8) = plot(T,log10(p_Tetot),'-k','DisplayName','FactSage Sum Te1-Te7');
for i = 1:8
    line(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({line(i).DisplayName},size(line(i).XData)));
end

TBrooks = T(T>800);
pTeTotBrooksmmHg = 10.^(-5960.2./TBrooks + 7.5999); 
pTeTotBrooksatm = pTeTotBrooksmmHg/760;
line(9) = plot(TBrooks,log10(pTeTotBrooksatm),'DisplayName','Brooks Total Te over Te');
line(9).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({line(9).DisplayName},size(line(9).XData)));

legend('Te_1', 'Te_2','Te_3','Te_4', 'Te_5','Te_6','Te_7','Te_{all}','Brooks Te_{all}')



% Te2 over CdTe
figure(2)
clf
hold on
ylim([-20 0])
title('Vapor Pressure of Te_2 over CdTe')
xlabel('T (K)')
ylabel('log_{10}(p in atm)')


% Te2 above CdTe!  
Tbrebrick = 1050:10:1230;
log10pTe2_Brebrick = (-10000./Tbrebrick + 6.346 );
line(10) = plot(Tbrebrick,log10pTe2_Brebrick,"DisplayName",'Brebrick Te2');
line(10).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({line(10).DisplayName},size(line(10).XData)));

% Te2 above CdTe!  
TLB1 = 300:10:1370;
TLB2 = 1370:10:1600;
pCdTeLB1 = 10.^(6.496 - 9580./TLB1);
pCdTeLB2 = 10.^(5.37 - 8050./TLB2);
pTe2LB1 = pCdTeLB1 / 3;
pTe2LB2 = pCdTeLB2 / 3;

line(12) = plot(TLB1,log10(pTe2LB1),"DisplayName",'LB T<Tm');
line(12).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({line(12).DisplayName},size(line(12).XData)));
line(13) = plot(TLB2,log10(pTe2LB2),"DisplayName",'LB T>Tm');
line(13).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series:',repmat({line(13).DisplayName},size(line(13).XData)));

legend('Brebrick 2010','LB T<T_m','LB T>T_m')

% plotting Gibbs energies per Te atom
figure(3)
clf
hold on
plot(T,G_Te_ls(T,1,1,P_units), T,G_Te_gv(T,1,1,P_units), T,0.5*G_Te2_gv(T,1,1,P_units) )
plot(T,(2*G_Te_ls(T,1,1,P_units)-G_Te2_gv(T,1,1,P_units))/2,'ko')
plot(T, 0.5*kB_eV*T.*log(p_Te2),'g-')
legend('G^o_{l,s}', 'G^o_{Te_vap}','G^o_{Te_2,vap}','\mu_{Te} from (2*G^o_{Te,ls} - G^o_{Te_2,vap})/2','\mu_{Te} from 0.5*kBT ln(p_{vap,Te_2})')
xlabel('T (K)')
ylabel('Gibbs Energy (eV/atom)')

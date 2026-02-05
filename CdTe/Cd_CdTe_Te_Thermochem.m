clear all
kB = 8.617333262e-5;

T = 300:1:1400;
T = T(:);


[G0_Cd_gv] = G0_Cd_gv(T, 1, 1, 'atm');
[G0_Cd_ls] = G0_Cd_ls(T, 1, 1, 'atm');
[G0_Te_gv] = G0_Te_gv(T, 1, 1, 'atm');
[G0_Te_ls] = G0_Te_ls(T, 1, 1, 'atm');
[G0_CdTe_ls] = G0_CdTe_ls(T, 1, 1, 'atm');

DG_Cd = G0_Cd_gv-G0_Cd_ls;
DG_Te = G0_Te_gv-G0_Te_ls;

pvapCd_over_Cd = 1*exp(-DG_Cd./(kB*T));
pvapTe2_over_Te = 1*exp(-DG_Te./(kB*T));

data_out = [T pvapCd_over_Cd pvapTe2_over_Te G0_Cd_ls G0_Cd_gv G0_Te_ls G0_Te_gv G0_CdTe_ls];



figure()
clf
hold on
plot(1./T,log10(pvapCd_over_Cd),'b-')
plot(1./T,log10(pvapTe2_over_Te),'r-')
ylim([-10 1])
legend(["pCd over Cd" "pTe2 over Te"])

figure()
clf
hold on
plot(T,G0_Cd_ls)
plot(T,G0_Cd_gv)
plot(T,G0_Te_ls)
plot(T,G0_Te_gv)
plot(T,G0_CdTe_ls)
legend(["Cd_{(s,l)}" "Cd_{(vap,1 atm)}" "Te_{(s,l)}" "Te_{(vap,1 atm)}" "CdTe_{(s,l)}"])
xlabel("T (K)")
ylabel("G_0 (eV/molecule)")

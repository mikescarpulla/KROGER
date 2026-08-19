for i=1:2
    figure(i)
    clf
    hold on
end

kB = 8.617E-05;

T_inv_Te_Terich = [0.73284
    0.73389
    0.73608
    0.73751
    0.73999
    0.74257
    0.74868
    0.755
    0.76
    0.76213
    0.764
    0.77
    0.78
    0.79
    0.8
    0.80994
    0.83
    0.85
    0.88083
    0.90078
    0.95
    1
    1.05
    1.08026
    1.20011
    1.219512195
    1.25
    1.282051282
    1.315789474
    1.351351351
    1.388888889];

log10_p_Te2_Terich = [-2.57643
    -2.23823
    -1.95153
    -1.85368
    -1.69746
    -1.59445
    -1.38501
    -1.22184875
    -1.124938737
    -1.07943
    -1.055517328
    -0.991399828
    -0.903089987
    -0.829738285
    -0.795880017
    -0.76011
    -0.725842151
    -0.725842151
    -0.81848
    -0.89745
    -1.086186148
    -1.337242168
    -1.585026652
    -1.78329
    -2.49231
    -2.601810335
    -2.765885476
    -2.934114856
    -3.106714315
    -3.283916965
    -3.465975081];

% linear subset
mask = T_inv_Te_Terich > 1.07;
T_low  = T_inv_Te_Terich(mask);
p_Te2_low = log10_p_Te2_Terich(mask);

% linear fit
p = polyfit(T_low, p_Te2_low, 1);

% extrapolate to T_m_Te = 722 K to get p_Te2 @T_melt
T_inv_target_722 = 1000/722;
log10_p_Te2_722 = polyval(p, T_inv_target_722)
p_Te2_722 = 10^log10_p_Te2_722

T_inv_Te_Terich = [T_inv_Te_Terich; T_inv_target_722];
log10_p_Te2_Terich = [log10_p_Te2_Terich; log10_p_Te2_722];
p_Te2_Terich = 10.^log10_p_Te2_Terich;


Te_poly_Terich = polyfit(T_inv_Te_Terich, log10_p_Te2_Terich, numel(T_inv_Te_Terich)-1); % use to calculate G_CdTe later on

figure(1)
plot(T_inv_Te_Terich, log10_p_Te2_Terich, 'DisplayName', 'pTe_2, Te-rich');

T_2_Te_Terich = 1000./T_inv_Te_Terich;
T_1_Te_Terich = 300:T_2_Te_Terich(end);

mu_Te_KROGER_Terich = G_Te_ls(T_1_Te_Terich, 1, 1, 'atm');
mu_Te_722 = mu_Te_KROGER_Terich(end);

% dmu_Te = dmu0_Te(722) + kBTln(P_Te/P_ref)/2


dmu_Te_Terich_with_pTe = zeros(numel(T_inv_Te_Terich),1);

for i = 1:numel(dmu_Te_Terich_with_pTe)
    dmu_Te_Terich_with_pTe(i) = mu_Te_722 + (kB * T_2_Te_Terich(i) * log(p_Te2_Terich(i)/p_Te2_722))/2 ;
end

T_Te_Terich = [T_1_Te_Terich'; T_2_Te_Terich];
dmu_Te_Terich = [mu_Te_KROGER_Terich; dmu_Te_Terich_with_pTe];

figure(2)
plot(T_Te_Terich,   dmu_Te_Terich,'or', 'DisplayName', '\Delta\mu_{Te}, Te-rich');







T_inv_Te_Cdrich = [0.73284
    0.73427
    0.73608
    0.73971
    0.74658
    0.75412
    0.7628
    0.7712
    0.78169
    0.78856
    0.80736];

log10_p_Te2_Cdrich = [-2.57643
    -2.90948
    -3.21506
    -3.6185
    -4.00992
    -4.37387
    -4.67945
    -4.96099
    -5.24597
    -5.43482
    -5.87259];

p_Te2_Cdrich = 10.^log10_p_Te2_Cdrich;
Te_poly_Cdrich = polyfit(T_inv_Te_Cdrich,log10_p_Te2_Cdrich,numel(T_inv_Te_Cdrich)-1);

figure(1)

title('log_{10}(p_{Te2}) vs 1000/T');
plot(T_inv_Te_Cdrich, log10_p_Te2_Cdrich, 'DisplayName', 'pTe_2, Cd-rich');
legend show;

T_Te_Cdrich = 1000./T_inv_Te_Cdrich;
mu_Te_KROGER_Cdrich = G_Te_ls(T_Te_Cdrich, 1, 1, 'atm');

% dmu_Te = dmu0_Te(722) + kBTln(P_Te/P_ref)/2

dmu_Te_Cdrich = zeros(numel(T_inv_Te_Cdrich),1);

for i = 1:numel(dmu_Te_Cdrich)
    dmu_Te_Cdrich(i) = mu_Te_722 + (kB * T_Te_Cdrich(i) * log(p_Te2_Cdrich(i)/p_Te2_722))/2 ;
end

figure(2)

title('dmu_{Te} vs T');
plot(T_Te_Cdrich, dmu_Te_Cdrich,'ob', 'DisplayName', '\Delta\mu_{Te}, Cd-rich');
legend show;

%% mu_Cd calc.
T_inv_Cd_Cdrich = [0.73241
    0.732600733
    0.733184315
    0.733917075
    0.734796926
    0.736268497
    0.742195197
    0.752237666
    0.759924207
    0.77166
    0.7829
    0.79966
    0.82017
    0.86614
    0.87992
    0.89447
    0.90869
    0.91916
    0.93085
    0.94463
    0.96051
    0.97396
    0.99148
    0.99942
    1.0993
    1.18077
    1.219512195
    1.25
    1.282051282
    1.315789474
    1.351351351
    1.388888889];

log10_p_Cd_Cdrich = [0.09057
    0.100370545
    0.235192448
    0.33157647
    0.407148462
    0.486809953
    0.629136055
    0.742901585
    0.80062406
    0.81588
    0.80962
    0.76105
    0.68429
    0.46497
    0.39918
    0.33338
    0.25662
    0.20336
    0.1501
    0.07334
    -0.00499
    -0.07549
    -0.16478
    -0.20081
    -0.73501
    -1.16268
    -1.339836835
    -1.480855665
    -1.625444969
    -1.773790256
    -1.92609188
    -2.082566666];

% linear subset
mask = T_inv_Cd_Cdrich > 1;
T_low  = T_inv_Cd_Cdrich(mask);
p_Cd_low = log10_p_Cd_Cdrich(mask);

% linear fit
p = polyfit(T_low, p_Cd_low, 1);

% extrapolate to T_m_Cd = 594 K
T_inv_target_594 = 1000/594;
log10_p_Cd_594 = polyval(p, T_inv_target_594)
p_Cd_594 = 10^log10_p_Cd_594
T_inv_Cd_Cdrich = [T_inv_Cd_Cdrich; T_inv_target_594];
log10_p_Cd_Cdrich = [log10_p_Cd_Cdrich; log10_p_Cd_594];
p_Cd_Cdrich = 10.^log10_p_Cd_Cdrich;

Cd_poly_Cdrich = polyfit(T_inv_Cd_Cdrich, log10_p_Cd_Cdrich, numel(T_inv_Cd_Cdrich)-1); % use to calculate G_CdTe later on

figure(1)
plot(T_inv_Cd_Cdrich, log10_p_Cd_Cdrich, 'DisplayName', 'pCd, Cd-rich');

T_2_Cd_Cdrich = 1000./T_inv_Cd_Cdrich;
T_1_Cd_Cdrich = 300:T_2_Cd_Cdrich(end);

mu_Cd_KROGER_Cdrich = G_Cd_ls(T_1_Cd_Cdrich, 1, 1, 'atm');
mu_Cd_594 = mu_Cd_KROGER_Cdrich(end);

% dmu_Te = dmu0_Te(722) + kBTln(P_Te/P_ref)/2

dmu_Cd_Cdrich_with_pTe = zeros(numel(T_inv_Cd_Cdrich),1);

for i = 1:numel(dmu_Cd_Cdrich_with_pTe)
    dmu_Cd_Cdrich_with_pTe(i) = mu_Cd_594 + (kB * T_2_Cd_Cdrich(i) * log(p_Cd_Cdrich(i)/p_Cd_594)) ;
end

T_Cd_Cdrich = [T_1_Cd_Cdrich' ;    T_2_Cd_Cdrich];

dmu_Cd_Cdrich = [mu_Cd_KROGER_Cdrich
    dmu_Cd_Cdrich_with_pTe];

figure(2)
plot(T_Cd_Cdrich, dmu_Cd_Cdrich,'or', 'DisplayName', '\Delta\mu_{Cd}, Cd-rich');

T_inv_Cd_Terich = [1.01783
    0.99942
    0.97848
    0.96095
    0.93713
    0.9217
    0.89954
    0.87948
    0.85445
    0.83395
    0.81432
    0.79503
    0.77838
    0.75865
    0.74785
    0.74035
    0.73903
    0.73825
    0.7366
    0.73274
    0.73241];

log10_p_Cd_Terich = [-4.71406
    -4.501
    -4.23312
    -4.01381
    -3.70519
    -3.52034
    -3.22896
    -2.95012
    -2.60547
    -2.31409
    -2.01802
    -1.71567
    -1.41802
    -1.00759
    -0.68174
    -0.41386
    -0.34337
    -0.29481
    -0.21178
    0.00127
    0.09057];

p_Cd_Terich = 10.^log10_p_Cd_Terich;


Cd_poly_Terich = polyfit(T_inv_Cd_Terich, log10_p_Cd_Terich, numel(T_inv_Cd_Terich)-1);

figure(1)
title('log_{10}(p_{Cd}) vs 1000/T')
plot(T_inv_Cd_Terich, log10_p_Cd_Terich, 'DisplayName', 'pCd, Te-rich');
legend show;

T_Cd_Terich = 1000./T_inv_Cd_Terich;

mu_Cd_KROGER_Terich = G_Cd_ls(T_Cd_Terich, 1, 1, 'atm');

% dmu_Te = dmu0_Te(722) + kBTln(P_Te/P_ref)/2

dmu_Cd_Terich = zeros(numel(T_inv_Cd_Terich),1);

for i = 1:numel(dmu_Cd_Terich)
    dmu_Cd_Terich(i) = mu_Cd_594 + (kB * T_Cd_Terich(i) * log(p_Cd_Terich(i)/p_Cd_594)) ;
end

figure(2)

title('dmu_{Cd} vs T')
plot(T_Cd_Terich, dmu_Cd_Terich,'ob', 'DisplayName', '\Delta\mu_{Cd}, Te-rich');
legend show;

%% all the plots before this part - just to make sure everything looks alright
%% fitting to get mu values from T inputs for further calculation in KROGER

% mu_Te at Te rich
[TT, idx] = sort(T_Te_Terich(:));
muT = dmu_Te_Terich(idx);
[TT, ia] = unique(TT, 'stable');      % important if any duplicate T exists
muT = muT(ia);
muTe_Terich = @(Tq) interp1(TT, muT, Tq, 'pchip');

% mu_Te at Cd rich
[TT, idx] = sort(T_Te_Cdrich(:));
muT = dmu_Te_Cdrich(idx);
[TT, ia] = unique(TT, 'stable');
muT = muT(ia);
muTe_Cdrich = @(Tq) interp1(TT, muT, Tq, 'pchip');

% mu_Cd at Cd rich
[TT, idx] = sort(T_Cd_Cdrich(:));
muT = dmu_Cd_Cdrich(idx);
[TT, ia] = unique(TT, 'stable');
muT = muT(ia);
muCd_Cdrich = @(Tq) interp1(TT, muT, Tq, 'pchip');

% mu_Cd at Te rich
[TT, idx] = sort(T_Cd_Terich(:));
muT = dmu_Cd_Terich(idx);
[TT, ia] = unique(TT, 'stable');
muT = muT(ia);
muCd_Terich = @(Tq) interp1(TT, muT, Tq, 'pchip');

figure(2)
% title('dmu_{Te} vs T (fitting)')
plot(305:10:1365, muTe_Terich(305:10:1365), 'DisplayName', '\Delta\mu_{Te}, Te-rich fitted')
plot(1238:1365, muTe_Cdrich(1238:1365), 'DisplayName', '\Delta\mu_{Te}, Cd-rich fitted')

% title('dmu_{Cd} vs T (fitting)')
plot(305:1365, muCd_Cdrich(305:1365), 'DisplayName', '\Delta\mu_{Cd}, Cd-rich fitted')
plot(982:1365, muCd_Terich(982:1365), 'DisplayName', '\Delta\mu_{Cd}, Te-rich fitted')
legend show;

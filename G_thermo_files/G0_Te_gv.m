function [G0_Te_gv] = G0_Te_gv(T, P_tot, X_i, P_units)

% define constants
q = 1.602176634e-19;
avo = 6.0221409e+23;
kB_eV = 8.617333262e-5;

% select the Pref for specified units of pressure 
if strcmp(P_units,'atm')
    P_ref = 1;
elseif strcmp(P_units,'Torr')
    P_ref = 760;   
elseif strcmp(P_units,'Bar')
    P_ref = 1;
elseif strcmp(P_units,'Pa')
    P_ref = 1e5;
else
    error('Units of pressure must be atm, Torr, Pa, or Bar, or you need to add more units with corresponding Pref in the file')
end

% work with T and P if they come in as col or row vectors
T = T(:);  % make T col vector
nT=size(T,1);
P_tot = P_tot(:)';  % make P a row vector
nP = size(P_tot,1);

% shape T and P to both be 2D arrays so that all calcs are vectorized (matrixized?)
T = T*ones(1,nP);   % copy T to make a matrix with as many columns as P has entries
P_tot = ones(nT,1)*P_tot;   % similarly for P - copy down as many rows as T has entries
G0_Te_gv = zeros(size(T));

% define masks needed based on Tranges given for the expressions for G0.  logical array same size as T with ones where true and 0 where false
mask1 = (T >= 298) & (T < 500);
mask2 = (T >= 500) & (T < 1000);
mask3 = (T >= 1000) & (T <= 3000);
G0_Te_gv = G0_Te_gv + mask1.*(203253.654 - 43.3817547*T - 20.7860000*T.*log(T));
G0_Te_gv = G0_Te_gv + mask2.*(202890.832 - 36.6298799*T + 1.455370000e-3*T.^2 + 22851.6000*T.^(-1) - 3.676316667e-7*T.^3 - 21.8727000*T.*log(T));
G0_Te_gv = G0_Te_gv + mask3.*(207495.815 - 82.4914885*T - 2.720890000e-3*T.^2 - 577900.000*T.^(-1) + 1.259381667e-7*T.^3 - 15.2801000*T.*log(T));

% now convert units to eV per Ga2O
G0_Te_gv = G0_Te_gv/(avo*q);   % eV/Ga2O molecule

% Now take Ptot and Xi into account.  
G0_Te_gv = G0_Te_gv + kB_eV*T.*( log(P_tot/P_ref) + log(X_i));

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen 
G0_Te_gv(G0_Te_gv==0) = Inf;
G0_Te_gv(isnan(G0_Te_gv)) = Inf;

end


% % View Data  Te     Units:  T(K) P(atm) Energy(J) Quantity(mol) 
% % Name: Tellurium
% % 
% %   G(T) J/mol - 1 atm  
% % 
% %              G(T)                  G(T)                   G(T)                     T(K)        
% % ____________ _____________________ ______________________ ________________________ ___________ 
% % 
% % S1         1 - 10544.6744          + 183.372880     T     + 1.583435900E-02 T^ 2   298 - 723   
% % S1         1 + 155014.500     T^-1 - 5.240416667E-06 T^ 3 - 35.6687000     T ln(T) 298 - 723   
% % S1         2 9160.59470            - 129.265375     T     - 3.623610000E-02 T^ 2   723 - 1150  
% % S1         2 - 1286810.00     T^-1 + 5.006366667E-06 T^ 3 + 13.0040000     T ln(T) 723 - 1150  
% % S1         3 - 12781.3493          + 174.901224     T     - 32.5596000     T ln(T) 1150 - 1600 
% % L1         4 - 17554.6235          + 685.877694     T     + 0.221943360     T^ 2   298 - 626   
% % L1         4 + 827927.650     T^-1 - 9.420746000E-05 T^ 3 - 126.318020     T ln(T) 298 - 626   
% % L1         5 - 3165763.16          + 46756.3503     T     + 7.09774900     T^ 2    626 - 723   
% % L1         5 + 258051100.     T^-1 - 1.306927800E-03 T^ 3 - 7196.40900     T ln(T) 626 - 723   
% % L1         6 180326.958            - 1500.57909     T     - 0.142016000     T^ 2   723 - 1150  
% % L1         6 - 24238450.0     T^-1 + 1.612973333E-05 T^ 3 + 202.743000     T ln(T) 723 - 1150  
% % L1         7 6328.68506            + 148.708305     T     - 32.5596000     T ln(T) 1150 - 1600 
% % G1         8 203253.654            - 43.3817547     T     - 20.7860000     T ln(T) 298 - 500   
% % G1         9 202890.832            - 36.6298799     T     + 1.455370000E-03 T^ 2   500 - 1000  
% % G1         9 + 22851.6000     T^-1 - 3.676316667E-07 T^ 3 - 21.8727000     T ln(T) 500 - 1000  
% % G1        10 207495.815            - 82.4914885     T     - 2.720890000E-03 T^ 2   1000 - 3000 
% % G1        10 - 577900.000     T^-1 + 1.259381667E-07 T^ 3 - 15.2801000     T ln(T) 1000 - 3000 
% % ____________ _____________________ ______________________ ________________________ ___________ 
% % 

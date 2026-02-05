function [G0_CdTe_ls] = G0_CdTe_ls(T, P_tot, X_i, P_units)

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
G0_CdTe_s= zeros(size(T));
G0_CdTe_l = zeros(size(T));
G0_CdTe_ls = zeros(size(T));


% define masks needed based on Tranges given for the expressions for G0.  logical array same size as T with ones where true and 0 where false

% solid
masks1 = (T>=298) .* (T<=723);      % solid
masks2 = (T>723) .* (T<=833);      % solid
masks3 = (T>833) .* (T<=1150);      % solid
masks4 = (T>1150) .* (T<=2000);      % solid
G0_CdTe_s = masks1.*(-114886.039 + 209.434630*T - 9.089000000E-3*T.^2 - 44.6350000*T.*log(T));
G0_CdTe_s = G0_CdTe_s + masks2.*(100923.919 - 2107.20451*T - 0.210370000*T.^2 - 24238450.0*T.^(-1) + 1.612973333E-5*T.^3 + 289.783000*T.*log(T));
G0_CdTe_s = G0_CdTe_s + masks3.*(59805.1601 - 1443.24743*T - 0.151105000*T.^2 - 24238450.0*T.^(-1) + 1.612973333E-5*T.^3 + 191.053000*T.*log(T));
G0_CdTe_s = G0_CdTe_s + masks4.*(-114193.112 + 206.039964*T - 9.089000000E-3*T.^2 - 44.2496000*T.*log(T));


% liquid
maskl1 = (T>298) .* (T<=400);      % liquid
maskl2 = (T>400) .* (T<=594);      % liquid
maskl3 = (T>594) .* (T<=626);      % liquid
maskl4 = (T>626) .* (T<=723);      % liquid
maskl5 = (T>723) .* (T<=1150);      % liquid
maskl6 = (T>1150) .* (T<=2000);      % liquid
G0_CdTe_l = maskl1.*(-93325.7245 + 780.722169*T + 0.215669452*T.^2 + 820961.290*T.^(-1) - 9.420746000E-5*T.^3 - 148.362261*T.*log(T));
G0_CdTe_l = G0_CdTe_l + maskl2.*(-70653.8162 + 320.466017*T + 0.106783443*T.^2 - 443887.796*T.^(-1) - 6.530767884E-5*T.^3 - 73.1866302*T.*log(T));
G0_CdTe_l = G0_CdTe_l + maskl3.*(-95623.0031 + 829.763994*T + 0.221943360*T.^2 + 827927.650*T.^(-1) - 9.420746000E-5*T.^3 - 156.024420*T.*log(T));
G0_CdTe_l = G0_CdTe_l + maskl4.*(-3243831.54 + 46900.2366*T + 7.09774900*T.^2 + 258051100.*T.^(-1) - 1.306927800E-3*T.^3 - 7226.11540*T.*log(T));
G0_CdTe_l = G0_CdTe_l + maskl5.*(102258.578 - 1356.69279*T - 0.142016000*T.^2 - 24238450.0*T.^(-1) + 1.612973333E-5*T.^3 + 173.036600*T.*log(T));
G0_CdTe_l = G0_CdTe_l + maskl6.*(-71739.6946 + 292.594604*T - 62.2660000*T.*log(T));


% choose the min Go between s and l to choose the right phase at each T
G0_CdTe_ls = min( cat(3, G0_CdTe_s, G0_CdTe_l),[],3);  % stack the two Go matrices along a dummy dimension, then take the min along that dimension.  Yes the [] is needed in the syntax for min() 

% now convert units to eV per Ga2O
G0_CdTe_ls = G0_CdTe_ls/(avo*q);   % eV/Ga2O molecule

% Now take Ptot and Xi into account.  For solids and liquids, it's only Xi
% that matters while gasses/vapors have the P/Pref term too.  
G0_CdTe_ls = G0_CdTe_ls + kB_eV*T.* log(X_i);

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen 
G0_CdTe_ls(G0_CdTe_ls==0) = Inf;

end



% % View Data  CdTe     Units:  T(K) P(atm) Energy(J) Quantity(mol) 
% % Name: Cadmium Telluride
% % 
% %   G(T) J/mol - 1 atm  
% % 
% %              G(T)                     G(T)                   G(T)                     T(K)        
% % ____________ ________________________ ______________________ ________________________ ___________ 
% % 
% % S1         1 - 114886.039             + 209.434630     T     - 9.089000000E-03 T^ 2   298 - 723   
% % S1         1 - 44.6350000     T ln(T)                                                 298 - 723   
% % S1         2 100923.919               - 2107.20451     T     - 0.210370000     T^ 2   723 - 833   
% % S1         2 - 24238450.0     T^-1    + 1.612973333E-05 T^ 3 + 289.783000     T ln(T) 723 - 833   
% % S1         3 59805.1601               - 1443.24743     T     - 0.151105000     T^ 2   833 - 1150  
% % S1         3 - 24238450.0     T^-1    + 1.612973333E-05 T^ 3 + 191.053000     T ln(T) 833 - 1150  
% % S1         4 - 114193.112             + 206.039964     T     - 9.089000000E-03 T^ 2   1150 - 2000 
% % S1         4 - 44.2496000     T ln(T)                                                 1150 - 2000 
% % L1         5 - 93325.7245             + 780.722169     T     + 0.215669452     T^ 2   298 - 400   
% % L1         5 + 820961.290     T^-1    - 9.420746000E-05 T^ 3 - 148.362261     T ln(T) 298 - 400   
% % L1         6 - 70653.8162             + 320.466017     T     + 0.106783443     T^ 2   400 - 594   
% % L1         6 - 443887.796     T^-1    - 6.530767884E-05 T^ 3 - 73.1866302     T ln(T) 400 - 594   
% % L1         7 - 95623.0031             + 829.763994     T     + 0.221943360     T^ 2   594 - 626   
% % L1         7 + 827927.650     T^-1    - 9.420746000E-05 T^ 3 - 156.024420     T ln(T) 594 - 626   
% % L1         8 - 3243831.54             + 46900.2366     T     + 7.09774900     T^ 2    626 - 723   
% % L1         8 + 258051100.     T^-1    - 1.306927800E-03 T^ 3 - 7226.11540     T ln(T) 626 - 723   
% % L1         9 102258.578               - 1356.69279     T     - 0.142016000     T^ 2   723 - 1150  
% % L1         9 - 24238450.0     T^-1    + 1.612973333E-05 T^ 3 + 173.036600     T ln(T) 723 - 1150  
% % L1        10 - 71739.6946             + 292.594604     T     - 62.2660000     T ln(T) 1150 - 1600 
% % ____________ ________________________ ______________________ ________________________ ___________ 
% % 

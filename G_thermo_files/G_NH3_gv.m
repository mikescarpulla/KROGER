function [G0_NH3_gv] = G0_NH3_gv(T, P_tot, X_i, P_units)
% gives the thermodynamic variables of NH3 vapor as function of T, P 
% T in K, P in atm or Torr (must put in string argument 'atm' or 'Torr' etc
%
% G computed in kJ/mol then after these are computed we
% convert to eV
% T and P should be vectors

% constants
q = 1.602176634e-19;
avo = 6.0221409e+23;
kB_eV = 8.617333262e-5;


if strcmp(P_units,'atm')
    P_ref = 1;
elseif strcmp(P_units,'Torr')
    P_ref = 760;   
elseif strcmp(P_units,'Bar')
    P_ref = 1;
elseif strcmp(P_units,'Pa')
    P_ref = 1e5;
else
    error('Units of pressure must be atm, Torr, Pa, or Bar unless the ln(P/Pref) term is changed')
end


T = T(:);  % make T col vector
nT=size(T,1);
P_tot = P_tot(:)';  % make P a row vector
nP = size(P_tot,1);

% note here T and P are vectors
T = T*ones(1,nP);   % copy T to make a matrix with as many columns as P has entries
P_tot = ones(nT,1)*P_tot;   % similarly for P - copy down as many rows as T has entries
G0_NH3_gv = zeros(size(T));

mask1 = (T>=298) .* (T<=1800);  % 298 to 1800     logical array same size as T with ones where true and 0 where false
mask2 = (T>1800) .* (T<=5100);   % 1800-5100
mask3 = (T>5100) .* (T<=6000);   % 5100-6000


% below Tmelt (zeros except for below Tmelt entires)
% T region 1
G0_NH3_gv = G0_NH3_gv + mask1.*( -40730.6224 + 17.0712882*T - 1.883012667E-02*(T.^2) - 269022.941*(T.^-1) + 1.254181808E-06*(T.^3) - 2879.47060*log(T) - 28.6984236*T.*log(T) );

% T region 2
G0_NH3_gv = G0_NH3_gv + mask2.*( -416631.426 + 1407.67938*T + 2.789060271E-03*(T.^2) + 13660528.3*(T.^-1) + 97923.2870*log(T) - 23552.1383*(T.^0.5) - 173.525464*T.*log(T) );

% T region 3
G0_NH3_gv = G0_NH3_gv + mask3.*( -113600.190 + 405.115613*T - 80.7510000*T.*log(T) );

% now convert units to eV per NH3 molecule
G0_NH3_gv = G0_NH3_gv/(avo*q);   % eV/NH3 molecule

% Now add in the pressure (zero for condensed phases) 
G0_NH3_gv = G0_NH3_gv + kB_eV*T.*( log(P_tot/P_ref) + log(X_i));

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen 
G0_NH3_gv(G0_NH3_gv==0) = Inf;


end

% View Data  NH3     Units:  T(K) P(atm) Energy(J) Quantity(mol) 
% 
% Name: Ammonia
% 
% 
% 
%   G(T) J/mol - 1 atm  
% 
% 
% 
%               G(T)                     G(T)                   G(T)                     T(K)        
% 
% _____________ ________________________ ______________________ ________________________ ___________ 
% 
% 
% 
% G1         1  - 40730.6224             + 17.0712882     T     - 1.883012667E-02 T^ 2   298 - 1800  
% 
% G1         1  - 269022.941     T^-1    + 1.254181808E-06 T^ 3 - 2879.47060     ln(T)   298 - 1800  
% 
% G1         1  - 28.6984236     T ln(T)                                                 298 - 1800  
% 
% G1         2  - 416631.426             + 1407.67938     T     + 2.789060271E-03 T^ 2   1800 - 5100 
% 
% G1         2  + 13660528.3     T^-1    + 97923.2870     ln(T) - 23552.1383     T^0.5   1800 - 5100 
% 
% G1         2  - 173.525464     T ln(T)                                                 1800 - 5100 
% 
% G1         3  - 113600.190             + 405.115613     T     - 80.7510000     T ln(T) 5100 - 6000 
% 
% Aq1         4 - 80291.0000             - 111.294000     T     Assume Cp = 0            298 - 298   
% 
% _____________ ________________________ ______________________ ________________________ ___________ 
% 
% 

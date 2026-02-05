function [G0_N2_gv] = G0_N2_gv(T, P_tot, X_i, P_units)
% gives the thermodynamic variables of N2 vapor as function of T, P  
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
G0_N2_gv = zeros(size(T));

mask1 = (T>=298) .* (T<=1600);  % 298 to 1600     logical array same size as T with ones where true and 0 where false
mask2 = (T>1600) .* (T<=6000);   % 1600-6000


% below Tmelt (zeros except for below Tmelt entires)
% T region 1
G0_N2_gv = G0_N2_gv + mask1.*( -20407.7348 - 80.6780615*T - 8.694894573E-03*(T.^2) + 115077.483*(T.^-1) + 7.529653737E-07*(T.^ 3) + 2939.80512*log(T) - 17.0812209*T.*log(T) );

% T region 2
G0_N2_gv = G0_N2_gv + mask2.*( -110841.938 + 147.166947*T + 5872430.64*(T.^-1) + 20259.1042*log(T) - 2702.05404*(T.^0.5) - 43.9005434*T.*log(T) );

% now convert units to eV per Ga atom
G0_N2_gv = G0_N2_gv/(avo*q);   % eV/Ga atom


% Now add in the pressure (zero for condensed phases) 
G0_N2_gv = G0_N2_gv + kB_eV*T.*( log(P_tot/P_ref) + log(X_i));

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen 
G0_N2_gv(G0_N2_gv==0) = Inf;


end




% View Data  N2     Units:  T(K) P(atm) Energy(J) Quantity(mol) 
% 
% Name: Nitrogen
% 
% 
% 
%   G(T) J/mol - 1 atm  
% 
% 
% 
%               G(T)                     G(T)                     G(T)                     T(K)        
% 
% _____________ ________________________ ________________________ ________________________ ___________ 
% 
% 
% 
% G1         1  - 20407.7348             - 80.6780615     T       - 8.694894573E-03 T^ 2   298 - 1600  
% 
% G1         1  + 115077.483     T^-1    + 7.529653737E-07 T^ 3   + 2939.80512     ln(T)   298 - 1600  
% 
% G1         1  - 17.0812209     T ln(T)                                                   298 - 1600  
% 
% G1         2  - 110841.938             + 147.166947     T       + 5872430.64     T^-1    1600 - 6000 
% 
% G1         2  + 20259.1042     ln(T)   - 2702.05404     T^0.5   - 43.9005434     T ln(T) 1600 - 6000 
% 
% Aq1         3 307132.388               - 4890.38302     T       - 0.758814500     T^ 2   298 - 523   
% 
% Aq1         3 - 24395850.0     T^-1    + 742.572000     T ln(T)                          298 - 523   
% 
% _____________ ________________________ ________________________ ________________________ ___________ 
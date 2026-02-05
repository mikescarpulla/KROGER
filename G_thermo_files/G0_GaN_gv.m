function [G0_GaN_gv] = G0_GaN_gv(T, P_tot, X_i, P_units)
% gives the thermodynamic variables of GaN as function of T, P 
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
nT = size(T,1);
P_tot = P_tot(:)';  % make P a row vector
nP = size(P_tot,1);


T = T*ones(1,nP);   % copy T to make a matrix with as many columns as P has entries
P_GaN = ones(nT,1)*P_GaN;   % similarly for P - copy down as many rows as T has entries
G0_GaN_gv = zeros(size(T));

mask1 = (T>298) .* (T<=1773);  % 0 to 1773     logical array same size as T with ones where true and 0 where false

% up to 900 K 
G0_GaN_gv = G0_GaN_gv + mask1.*(175728.000  - 225.936000.*T); 

% now convert units to eV per GaN
G0_GaN_gv = G0_GaN_gv/(avo*q);   % eV/GaN molecule

% Now take Ptot and Xi into account.  For solids and liquids, it's only Xi
% that matters while gasses/vapors have the P/Pref term too.
G0_GaN_gv = G0_GaN_gv + kB_eV*T.*(log(P_tot/P_ref) + log(X_GaN));

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen 
G0_GaN_gv(G0_GaN_gv==0) = Inf;

end


% % View Data  GaN     Units:  T(K) P(atm) Energy(J) Quantity(mol) 
% % Name: Gallium Nitride
% % 
% %   G(T) J/mol - 1 atm  
% % 
% %              G(T)                     G(T)               G(T)                   T(K)       
% % ____________ ________________________ __________________ ______________________ __________ 
% % 
% % S1         1 - 121372.707             + 227.983015     T - 4.497800000E-03 T^ 2 298 - 1773 
% % S1         1 - 38.0744000     T ln(T)                                           298 - 1773 
% % G1         2 175728.000               - 225.936000     T                        298 - 298  
% % ____________ ________________________ __________________ ______________________ __________ 
% % 

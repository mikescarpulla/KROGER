function [G0_Ar_gv] = G0_Ar_gv(T, P_tot, X_i, P_units)
% gives the thermodynamic variables of Ar

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
G0_Ar_gv = zeros(size(T));

mask1 = (T>298) .* (T<=6000);  % 0 to 600     logical array same size as T with ones where true and 0 where false


% T region 1
G0_Ar_gv = G0_Ar_gv + mask1.*( -6197.34590 - 15.5193109*T - 20.7860000*T.*log(T) );


% now convert units to eV per Ga atom
G0_Ar_gv = G0_Ar_gv/(avo*q);   % eV/Ga atom


% Now add in the pressure (zero for condensed phases) 
G0_Ar_gv = G0_Ar_gv + kB_eV*T.*( log(P_tot/P_ref) + log(X_i));

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen 
G0_Ar_gv(G0_Ar_gv==0) = Inf;


end



% % % 
% % % 
% % % View Data  Ar     Units:  T(K) P(atm) Energy(J) Quantity(mol) 
% % % Name: Argon
% % % 
% % %   G(T) J/mol - 1 atm  
% % % 
% % %               G(T)         G(T)               G(T)                     T(K)       
% % % _____________ ____________ __________________ ________________________ __________ 
% % % 
% % % G1         1  - 6197.34590 - 15.5193109     T - 20.7860000     T ln(T) 298 - 6000 
% % % Aq1         2 - 12134.0000 - 59.4130000     T Assume Cp = 0            298 - 298  
% % % _____________ ____________ __________________ ________________________ __________ 
% % % 

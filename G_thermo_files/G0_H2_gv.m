function [G0_H2_gv] = G0_H2_gv(T, P_tot, X_i, P_units)
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
G0_H2_gv = zeros(size(T));

mask1 = (T>=298) .* (T<=6000);  % 0 to 6000     logical array same size as T with ones where true and 0 where false

% T region 1
G0_H2_gv = G0_H2_gv + mask1.*(211801.654 + 24.6096891*T - 20.7860000*T.*log(T));

% now convert units to eV per H2 molecule
G0_H2_gv = G0_H2_gv/(avo*q);   % eV/H2 molecule


% Now add in the pressure (zero for condensed phases) 
G0_H2_gv = G0_H2_gv + kB_eV*T.*( log(P_tot/P_ref) + log(X_i));

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen 
G0_H2_gv(G0_H2_gv==0) = Inf;


end




% View Data  H     Units:  T(K) P(atm) Energy(J) Quantity(mol) 
% 
% Name: Hydrogen
% 
% 
% 
%   G(T) J/mol - 1 atm  
% 
% 
% 
%              G(T)       G(T)               G(T)                     T(K)       
% 
% ____________ __________ __________________ ________________________ __________ 
% 
% 
% 
% G1         1 211801.654 + 24.6096891     T - 20.7860000     T ln(T) 298 - 6000 
% 
% ____________ __________ __________________ ________________________ __________
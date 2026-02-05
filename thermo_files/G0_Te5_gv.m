function [G0_Te5_gv] = G0_Te5_gv(T, P_tot, X_i, P_units)

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
G0_Te5_gv = zeros(size(T));

% define masks needed based on Tranges given for the expressions for G0.  logical array same size as T with ones where true and 0 where false
mask1 = (T >= 298) & (T <= 1500);
G0_Te5_gv = G0_Te5_gv + mask1 .* (187116.902 + 262.128648*T - 1.964374304e-4*T.^2 + 107727.386*T.^(-1) + 5.160401344e-8*T.^3 - 6.972537041e-12*T.^4 - 107.906525*T.*log(T));


% now convert units to eV per Ga2O
G0_Te5_gv = G0_Te5_gv/(avo*q);   % eV/Ga2O molecule

% Now take Ptot and Xi into account.
G0_Te5_gv = G0_Te5_gv + kB_eV*T.*( log(P_tot/P_ref) + log(X_i));

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen
G0_Te5_gv(G0_Te5_gv==0) = Inf;
G0_Te5_gv(isnan(G0_Te5_gv)) = Inf;

end


% % View Data  Te5     Units:  T(K) P(atm) Energy(J) Quantity(mol) 
% % Name: tellurium
% % 
% %   G(T) J/mol - 1 atm  
% % 
% %              G(T)                     G(T)                   G(T)                   T(K)       
% % ____________ ________________________ ______________________ ______________________ __________ 
% % 
% % G1         1 187116.902               + 262.128648     T     - 1.964374304E-04 T^ 2 298 - 1500 
% % G1         1 + 107727.386     T^-1    + 5.160401344E-08 T^ 3 - 6.972537041E-12 T^ 4 298 - 1500 
% % G1         1 - 107.906525     T ln(T)                                               298 - 1500 
% % ____________ ________________________ ______________________ ______________________ __________ 
% % 

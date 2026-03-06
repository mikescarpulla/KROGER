function [G0_As3_gv] = G0_As3_gv(T, P_tot, X_i, P_units)
% T in kelvin
% Ptot in P_units - 'atm' for example 
% X_i = mole fraction, which for gases is partial pressure/total pressure

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
G0_As3_gv= zeros(size(T));


% define masks needed based on Tranges given for the expressions for G0.  logical array same size as T with ones where true and 0 where false

% gas1
mask = (T>298) .* (T<=1200);
G0_As3_gv = mask.*(241969.032 + 107.271661* T - 1.004160000e-04* T.^(2) + 138490.400* T.^(-1) - 62.0947440*T.* log(T));

% now convert units to eV per Cd
G0_As3_gv = G0_As3_gv/(avo*q);   % eV/Cd

% Now take Ptot and Xi into account.  For solids and liquids, it's only Xi
% that matters while gasses/vapors have the P/Pref term too.  
G0_As3_gv = G0_As3_gv + kB_eV*T.* ( log(P_tot/P_ref) + log(X_i));

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen 
G0_As3_gv(G0_As3_gv==0) = Inf;
G0_As3_gv(isnan(G0_As3_gv)) = Inf;

end



%%
% View Data   As3     Units:  T(K) P(atm) Energy(J) Quantity(mol) 
% Name: As3
% 
%   G(T) J/mol - 1 atm  
% 
%              G(T)                  G(T)                     G(T)                   T(K)       
% ____________ _____________________ ________________________ ______________________ __________ 
% 
% G1         1 241969.032            + 107.271661     T       - 1.004160000E-04 T^ 2 298 - 1200 
% G1         1 + 138490.400     T^-1 - 62.0947440     T ln(T)                        298 - 1200 
% ____________ _____________________ ________________________ ______________________ __________ 
% 

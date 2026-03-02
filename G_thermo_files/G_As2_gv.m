function [G0_As2_gv] = G0_As2_gv(T, P_tot, X_i, P_units)
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
G0_As2_gv= zeros(size(T));


% define masks needed based on Tranges given for the expressions for G0.  logical array same size as T with ones where true and 0 where false

% gas1
mask = (T>298) .* (T<=1200);
G0_As2_gv = mask.*(209223.218 + 9.44990604* T- 7.531200000e-05* T.^(2) + 101043.600 *T.^(-1) - 37.1999440* T.* log(T));

% now convert units to eV per Cd
G0_As2_gv = G0_As2_gv/(avo*q);   % eV/Cd

% Now take Ptot and Xi into account.  For solids and liquids, it's only Xi
% that matters while gasses/vapors have the P/Pref term too.  
G0_As2_gv = G0_As2_gv + kB_eV*T.* ( log(P_tot/P_ref) + log(X_i));

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen 
G0_As2_gv(G0_As2_gv==0) = Inf;
G0_As2_gv(isnan(G0_As2_gv)) = Inf;

end



%%
% View Data   As2     Units:  T(K) P(atm) Energy(J) Quantity(mol) 
% Name: As2
% 
%   G(T) J/mol - 1 atm  
% 
%              G(T)                  G(T)                     G(T)                   T(K)       
% ____________ _____________________ ________________________ ______________________ __________ 
% 
% G1         1 209223.218            + 9.44990604     T       - 7.531200000E-05 T^ 2 298 - 1200 
% G1         1 + 101043.600     T^-1 - 37.1999440     T ln(T)                        298 - 1200 
% ____________ _____________________ ________________________ ______________________ __________ 
% 

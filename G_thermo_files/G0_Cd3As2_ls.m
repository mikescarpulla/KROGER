function [G0_Cd3As2] = G0_Cd3As2_ls(T, P_tot, X_i, P_units)

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
G0_Cd3As2= zeros(size(T));


% define masks needed based on Tranges given for the expressions for G0.  logical array same size as T with ones where true and 0 where false

% solid
mask1 = (T>=298) .* (T<=994);      % solid
G0_Cd3As2= mask1.*(- 87283.0042 + 716.088475* T  - 5.962200000e-03*T.^ 2+ 642244.000 * T.^(-1) - 136.189200*T.* log(T));


% % choose the min Go between s and l to choose the right phase at each T
% G0_Cd3As2 = G0_Cd3As2_s;  % stack the two Go matrices along a dummy dimension, then take the min along that dimension.  Yes the [] is needed in the syntax for min() 

% now convert units to eV per Ga2O
G0_Cd3As2 = G0_Cd3As2/(avo*q);   % eV / FU

% Now take Ptot and Xi into account.  For solids and liquids, it's only Xi
% that matters while gasses/vapors have the P/Pref term too.  
G0_Cd3As2 = G0_Cd3As2 + kB_eV*T.* log(X_i);

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen 
G0_Cd3As2(G0_Cd3As2==0) = Inf;

end


% % 
% % View Data   Cd3As2     Units:  T(K) P(atm) Energy(J) Quantity(mol) 
% % Name: Cadmium Arsenide
% % 
% %   G(T) J/mol - 1 atm  
% % 
% %              G(T)                  G(T)                     G(T)                   T(K)      
% % ____________ _____________________ ________________________ ______________________ _________ 
% % 
% % S1         1 - 87283.0042          + 716.088475     T       - 5.962200000E-03 T^ 2 298 - 994 
% % S1         1 + 642244.000     T^-1 - 136.189200     T ln(T)                        298 - 994 
% % ____________ _____________________ ________________________ ______________________ _________ 
% % 

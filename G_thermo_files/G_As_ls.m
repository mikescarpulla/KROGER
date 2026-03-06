function [G0_As_ls] = G0_As_ls(T, P_tot, X_i, P_units)

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
G0_As_s= zeros(size(T));
G0_As_l = zeros(size(T));
G0_As_ls = zeros(size(T));


% define masks needed based on Tranges given for the expressions for G0.  logical array same size as T with ones where true and 0 where false

% solid
masks1 = (T>=298) .* (T<=1200);      % solid
G0_As_s = masks1.*(- 7270.44945  + 122.211126* T - 2.716149246e-03* T.^ 2 + 11599.5574* T.^(-1)- 23.3144103* T.* log(T)); 

% liquid
maskl1 = (T>298) .* (T<=500);      % liquid
maskl2 = (T>500) .* (T<=1090);      % liquid
maskl3 = (T>1090) .* (T<=1400);      % liquid
G0_As_l = maskl1.*(8743.13873 + 109.124595 *T - 2.716149246e-03* T.^(2) + 11599.5574* T.^(-1) - 23.3144103* T.* log(T));
G0_As_l = G0_As_l + maskl2.*(161061.021 - 2475.72395* T - 0.510986830* T.^(2) - 10663549.5* T.^(-1) + 1.437394233e-04* T.^(3) - 2.083333333e-08* T.^(4) + 385.998291 *T.* log(T));
G0_As_l = G0_As_l + maskl3.*(- 5259.40172 + 269.082850* T - 45.0000000* T.* log(T));

% choose the min Go between s and l to choose the right phase at each T
G0_As_ls = min( cat(3, G0_As_s, G0_As_l),[],3);  % stack the two Go matrices along a dummy dimension, then take the min along that dimension.  Yes the [] is needed in the syntax for min() 

% now convert units to eV per Ga2O
G0_As_ls = G0_As_ls/(avo*q);   % eV/Ga2O molecule

% Now take Ptot and Xi into account.  For solids and liquids, it's only Xi
% that matters while gasses/vapors have the P/Pref term too.  
G0_As_ls = G0_As_ls + kB_eV*T.* log(X_i);

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen 
G0_As_ls(G0_As_ls==0) = Inf;

end



%%
% View Data   As     Units:  T(K) P(atm) Energy(J) Quantity(mol) 
% Name: Arsenic
% 
%   G(T) J/mol - 1 atm  
% 
%              G(T)                     G(T)                     G(T)                     T(K)        
% ____________ ________________________ ________________________ ________________________ ___________ 
% 
% S1         1 - 7270.44945             + 122.211126     T       - 2.716149246E-03 T^ 2   298 - 1200  
% S1         1 + 11599.5574     T^-1    - 23.3144103     T ln(T)                          298 - 1200  
% L1         2 8743.13873               + 109.124595     T       - 2.716149246E-03 T^ 2   298 - 500   
% L1         2 + 11599.5574     T^-1    - 23.3144103     T ln(T)                          298 - 500   
% L1         3 161061.021               - 2475.72395     T       - 0.510986830     T^ 2   500 - 1090  
% L1         3 - 10663549.5     T^-1    + 1.437394233E-04 T^ 3   - 2.083333333E-08 T^ 4   500 - 1090  
% L1         3 + 385.998291     T ln(T)                                                   500 - 1090  
% L1         4 - 5259.40172             + 269.082850     T       - 45.0000000     T ln(T) 1090 - 1400 
% G1         5 295552.701               - 34.9928723     T       - 20.7861120     T ln(T) 298 - 1200  
% ____________ ________________________ ________________________ ________________________ ___________ 
% 

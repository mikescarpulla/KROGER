function [G0_As2Te3_ls] = G0_As2Te3_ls(T, P_tot, X_i, P_units)

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
G0_As2Te3_s= zeros(size(T));
G0_As2Te3_l = zeros(size(T));
G0_As2Te3_ls = zeros(size(T));


% define masks needed based on Tranges given for the expressions for G0.  logical array same size as T with ones where true and 0 where false

% solid
masks1 = (T>=298) .* (T<=750);      % solid
G0_As2Te3_s = masks1.*(- 86169.0053 + 702.742372 * T - 2.217520000e-02* T.^(2) + 929684.800* T.^(-1) - 135.185040* T.* log(T) ); 

% liquid
maskl1 = (T>298) .* (T<=500);      % liquid
maskl2 = (T>500) .* (T<=800);      % liquid
G0_As2Te3_l = maskl1.*(- 41151.5867 + 633.709261* T - 2.217520000e-02 *T.^(2) + 929685.000 *T.^(-1) - 135.185040* T.* log(T));
G0_As2Te3_l = G0_As2Te3_l + maskl2.*(- 47976.5267 + 839.945048* T - 167.360000* T.* log(T));

% choose the min Go between s and l to choose the right phase at each T
G0_As2Te3_ls = min( cat(3, G0_As2Te3_s, G0_As2Te3_l),[],3);  % stack the two Go matrices along a dummy dimension, then take the min along that dimension.  Yes the [] is needed in the syntax for min() 

% now convert units to eV per Ga2O
G0_As2Te3_ls = G0_As2Te3_ls/(avo*q);   % eV/Ga2O molecule

% Now take Ptot and Xi into account.  For solids and liquids, it's only Xi
% that matters while gasses/vapors have the P/Pref term too.  
G0_As2Te3_ls = G0_As2Te3_ls + kB_eV*T.* log(X_i);

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen 
G0_As2Te3_ls(G0_As2Te3_ls==0) = Inf;

end






%%
% View Data   As2Te3     Units:  T(K) P(atm) Energy(J) Quantity(mol) 
% Name: Arsenic Tritelluride
% 
%   G(T) J/mol - 1 atm  
% 
%              G(T)                  G(T)                     G(T)                     T(K)      
% ____________ _____________________ ________________________ ________________________ _________ 
% 
% S1         1 - 86169.0053          + 702.742372     T       - 2.217520000E-02 T^ 2   298 - 750 
% S1         1 + 929684.800     T^-1 - 135.185040     T ln(T)                          298 - 750 
% L1         2 - 41151.5867          + 633.709261     T       - 2.217520000E-02 T^ 2   298 - 500 
% L1         2 + 929685.000     T^-1 - 135.185040     T ln(T)                          298 - 500 
% L1         3 - 47976.5267          + 839.945048     T       - 167.360000     T ln(T) 500 - 800 
% ____________ _____________________ ________________________ ________________________ _________ 
% 

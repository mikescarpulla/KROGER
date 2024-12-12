function [G0_Te2_gv] = G0_Te2_gv(T, P_tot, X_i, P_units)

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
G0_Te2_gv = zeros(size(T));

% define masks needed based on Tranges given for the expressions for G0.  logical array same size as T with ones where true and 0 where false
mask1 = (T >= 298) & (T < 500);
mask2 = (T >= 500) & (T < 1000);
mask3 = (T >= 1000) & (T < 2100);
mask4 = (T >= 2100) & (T <= 3000);
G0_Te2_gv = G0_Te2_gv + mask1 .* (151341.335 - 19.4485261*T + 7.535000000e-4*T.^2 - 2.725716667e-6*T.^3 - 35.7007000*T.*log(T)); 
G0_Te2_gv = G0_Te2_gv + mask2 .* (155262.717 - 103.210571*T - 2.018790000e-2*T.^2 - 200257.000*T.^(-1) + 3.106300000e-6*T.^3 - 21.9053000*T.*log(T));
G0_Te2_gv = G0_Te2_gv + mask3 .* (133634.827 + 144.999946*T + 7.022800000e-3*T.^2 + 2176970.00*T.^(-1) - 6.666183333e-7*T.^3 - 58.4436000*T.*log(T));
G0_Te2_gv = G0_Te2_gv + mask4 .* (152003.363 - 3.80338102*T - 1.801695000e-3*T.^2 - 38.0321000*T.*log(T));

% now convert units to eV per Ga2O
G0_Te2_gv = G0_Te2_gv/(avo*q);   % eV/Ga2O molecule

% Now take Ptot and Xi into account.
G0_Te2_gv = G0_Te2_gv + kB_eV*T.*( log(P_tot/P_ref) + log(X_i));

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen
G0_Te2_gv(G0_Te2_gv==0) = Inf;
G0_Te2_gv(isnan(G0_Te2_gv)) = Inf;

end

% % View Data  Te2     Units:  T(K) P(atm) Energy(J) Quantity(mol) 
% % Name: Te2
% % 
% %   G(T) J/mol - 1 atm  
% % 
% %              G(T)                     G(T)                     G(T)                     T(K)        
% % ____________ ________________________ ________________________ ________________________ ___________ 
% % 
% % G1         1 151341.335               - 19.4485261     T       + 7.535000000E-04 T^ 2   298 - 500   
% % G1         1 - 2.725716667E-06 T^ 3   - 35.7007000     T ln(T)                          298 - 500   
% % G1         2 155262.717               - 103.210571     T       - 2.018790000E-02 T^ 2   500 - 1000  
% % G1         2 - 200257.000     T^-1    + 3.106300000E-06 T^ 3   - 21.9053000     T ln(T) 500 - 1000  
% % G1         3 133634.827               + 144.999946     T       + 7.022800000E-03 T^ 2   1000 - 2100 
% % G1         3 + 2176970.00     T^-1    - 6.666183333E-07 T^ 3   - 58.4436000     T ln(T) 1000 - 2100 
% % G1         4 152003.363               - 3.80338102     T       - 1.801695000E-03 T^ 2   2100 - 3000 
% % G1         4 - 38.0321000     T ln(T)                                                   2100 - 3000 
% % ____________ ________________________ ________________________ ________________________ ___________ 
% % 
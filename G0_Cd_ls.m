function [G0_Cd_ls] = G0_Cd_ls(T, P_tot, X_i, P_units)

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
G0_Cd_s= zeros(size(T));
G0_Cd_l = zeros(size(T));
G0_Cd_ls = zeros(size(T));


% define masks needed based on Tranges given for the expressions for G0.  logical array same size as T with ones where true and 0 where false

% solid
masks1 = (T>=298) .* (T<=594);      % solid
masks2 = (T>594) .* (T<=1500);      % solid
masks3 = (T>1500) .* (T<=1600);      % solid
G0_Cd_s = masks1.*(-7083.46898 + 99.5061986*T - 6.273908000e-3* T.^2 - 6966.36*T.^(-1) - 22.0442408*T.*log(T)); 
G0_Cd_s = G0_Cd_s + masks2.*(-20064.9716 + 256.812234*T + 8.832011496e-3*T.^2 + 1241289.61*T.^(-1) - 8.996036074e-7*T.^3 - 45.1611543*T.*log(T));
G0_Cd_s = G0_Cd_s + masks3.*(-9027.48876 + 148.205481*T - 29.7064000*T.*log(T));

% liquid
maskl1 = (T>298) .* (T<=400);      % liquid
maskl2 = (T>400) .* (T<=594);      % liquid
maskl3 = (T>594) .* (T<=1600);      % liquid
G0_Cd_l = maskl1.*(-955.024687 + 89.2092822*T - 6.273908000e-3*T.^2 - 6966.36000*T.^(-1) - 22.0442408*T.*log(T));
G0_Cd_l = G0_Cd_l + maskl2.*(21716.8836 - 371.046869*T - 0.115159917*T.^2 - 1271815.45*T.^(-1) + 2.889978116e-5*T.^3 + 53.1313898*T.*log(T));
G0_Cd_l = G0_Cd_l + maskl3.*(-3252.30331 + 138.251107*T - 29.7064000*T.*log(T));

% choose the min Go between s and l to choose the right phase at each T
G0_Cd_ls = min( cat(3, G0_Cd_s, G0_Cd_l),[],3);  % stack the two Go matrices along a dummy dimension, then take the min along that dimension.  Yes the [] is needed in the syntax for min() 

% now convert units to eV per Ga2O
G0_Cd_ls = G0_Cd_ls/(avo*q);   % eV/Ga2O molecule

% Now take Ptot and Xi into account.  For solids and liquids, it's only Xi
% that matters while gasses/vapors have the P/Pref term too.  
G0_Cd_ls = G0_Cd_ls + kB_eV*T.* log(X_i);

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen 
G0_Cd_ls(G0_Cd_ls==0) = Inf;

end



% % View Data  Cd     Units:  T(K) P(atm) Energy(J) Quantity(mol) 
% % Name: Cadmium
% % 
% %   G(T) J/mol - 1 atm  
% % 
% %              G(T)                  G(T)                     G(T)                     T(K)        
% % ____________ _____________________ ________________________ ________________________ ___________ 
% % 
% % S1         1 - 7083.46898          + 99.5061986     T       - 6.273908000E-03 T^ 2   298 - 594   
% % S1         1 - 6966.36000     T^-1 - 22.0442408     T ln(T)                          298 - 594   
% % S1         2 - 20064.9716          + 256.812234     T       + 8.832011496E-03 T^ 2   594 - 1500  
% % S1         2 + 1241289.61     T^-1 - 8.996036074E-07 T^ 3   - 45.1611543     T ln(T) 594 - 1500  
% % S1         3 - 9027.48876          + 148.205481     T       - 29.7064000     T ln(T) 1500 - 1600 
% % L1         4 - 955.024687          + 89.2092822     T       - 6.273908000E-03 T^ 2   298 - 400   
% % L1         4 - 6966.36000     T^-1 - 22.0442408     T ln(T)                          298 - 400   
% % L1         5 21716.8836            - 371.046869     T       - 0.115159917     T^ 2   400 - 594   
% % L1         5 - 1271815.45     T^-1 + 2.889978116E-05 T^ 3   + 53.1313898     T ln(T) 400 - 594   
% % L1         6 - 3252.30331          + 138.251107     T       - 29.7064000     T ln(T) 594 - 1600  
% % G1         7 105599.101            - 28.4149645     T       - 20.7861120     T ln(T) 298 - 1500  
% % ____________ _____________________ ________________________ ________________________ ___________ 
% % 

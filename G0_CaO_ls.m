function [G0_CaO_ls] = G0_CaO_ls(T, P_tot, X_i, P_units)
%
%  each substance should have a function Go_substance_ls.m for condensed
%  phases, and Go_substance_gv.m for gas/vapor.  See FAQs below.  
% 
% this function gives the standard Gibbs energy Go of substance as function of T, Ptot, and mole fraction which for a gas or ideal solution is X_i=P_partial/P_tot
% T and Ptot can come in as vectors, so output in general will be a 2D array.  
% T in K, P in atm, Torr etc (must put in string argument 'atm' or 'Torr' etc)
% G computed in kJ/mol becasue tables are this way then we convert to eV/formula unit at the end
% Logical masks are used like step functions vs T, allowing the whole expression to be built up piecewise over the list of T values  
%
% for a gas/vapor or mechanical mixture of phases such as ideal solution, mu(T,Ptot,Xi) = {Go(T,Pref) from Showmate equation} + kBT[ln(Ptot/Pref) + ln(X_i)]
% the total Gibbs energy of a mixture is Gmixture = Sum(Xi*mu_i) . Note how you end up with the xlnx terms this way, giving the ideal entropy of mixing
%
% FAQ1: Yes, we need to keep the P dependence in the input arguments for condensed phases becasue it makes the output an array of the right size even if the columns are identical 
% 
% FAQ2: why not put all phases in 1 file (l, s, and gas)?  Splitting out the gas/vapor phases can allow us to use these same files to compute the equilibrium vapor pressure for specifying thermodynamic conditions and things like that, rather than just having
% these files only spit out mu values.  
%
% FAQ3: waht if there are multiple allotropes for condensed phases?  List their Showmate equations (including any P dependencies).  Calculate the Go for all
% of them.  Then at the end of the section, just pick the minimum value for each T,P_tot point - that automatically gives us the Go for the stable phase under those conditions.  
% A future improvement could be to also spit out an array telling which phase is stable where.  Could be done with T,P masking and a unique ID for each phase.   
%
% FAQ4: why keep the  
%
% Handy preformated Showmate templates
% Go = mask1.*(A + B*T.^(-2) + C*T.^(-1) + D*T.^(0.5) + E*T + F*T.^(2) + G*T.^(3) + H*T.^(4) + J*T.*log(T));
% Go = Go + mask2.*(A + B*T.^(-2) + C*T.^(-1) + D*T.^(0.5) + E*T + F*T.^(2) + G*T.^(3) + H*T.^(4) + J*T.*log(T));
% Go = Go + mask3.*(A + B*T.^(-2) + C*T.^(-1) + D*T.^(0.5) + E*T + F*T.^(2) + G*T.^(3) + H*T.^(4) + J*T.*log(T));
% Go = Go + mask4.*(A + B*T.^(-2) + C*T.^(-1) + D*T.^(0.5) + E*T + F*T.^(2) + G*T.^(3) + H*T.^(4) + J*T.*log(T));

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
G0_CaO_ls = zeros(size(T));

% solid
mask1 = (T>0) .* (T<=2845);      %logical array same size as T with ones where true and 0 where false
G0_CaO_ls = G0_CaO_ls + mask1.*(- 651262.658 + 573572.991*T.^(-1) - 17163131.3*T.^(-2) - 535.615998*T.^(0.5) + 376.676564*T - 58.7911706*T.*log(T));

% liquid
mask2 = (T>298) .* (T<=2845);      % logical array same size as T with ones where true and 0 where false
G0_CaO_ls = G0_CaO_ls + mask2.*(- 571766.658 - 17163131.3*T.^(-2) + 573572.991*T.^(-1) + 348.735802*T - 535.615998*T.^(0.5) - 58.7911706*T.*log(T));



% now convert units to eV per Ga2O
G0_CaO_ls = G0_CaO_ls/(avo*q);   % eV/Ga2O molecule

% Now take Ptot and Xi into account.  
G0_CaO_ls = G0_CaO_ls + kB_eV*T.* log(X_i);

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen 
G0_CaO_ls(G0_CaO_ls==0) = Inf;

end


% % % 
% % % View Data  Fe3O4     Units:  T(K) P(atm) Energy(J) Quantity(mol) 
% % % Name: Ferrous ferric oxide
% % % 
% % %   G(T) J/mol - 1 atm  
% % % 
% % %              G(T)                  G(T)                     G(T)                     T(K)        
% % % ____________ _____________________ ________________________ ________________________ ___________ 
% % % 
% % % S1         1 - 1200277.81          + 1282.95850     T       - 1.901728207E-02 T^ 2   298 - 848   
% % % S1         1 + 3621668.52     T^-1 + 3.967903931E-05 T^ 3   - 3.104596095E-08 T^ 4   298 - 848   
% % % S1         1 - 110726059.     T^-2 - 207.930830     T ln(T)                          298 - 848   
% % % S1         2 - 1186819.15          + 1260.56231     T       + 3621668.52     T^-1    848 - 1870  
% % % S1         2 - 110726059.     T^-2 - 207.930830     T ln(T)                          848 - 1870  
% % % S2         3 - 1185253.60          + 1258.71640     T       + 3621668.52     T^-1    848 - 1870  
% % % S2         3 - 110726059.     T^-2 - 207.930830     T ln(T)                          848 - 1870  
% % % S2         4 - 1191672.58          + 1304.25027     T       - 213.384000     T ln(T) 1870 - 2500 
% % % S3         5 - 1030731.16          + 1253.59865     T       - 1.901728207E-02 T^ 2   298 - 848   
% % % S3         5 + 3621668.31     T^-1 + 3.967903931E-05 T^ 3   - 3.104596095E-08 T^ 4   298 - 848   
% % % S3         5 - 110726059.     T^-2 - 207.930830     T ln(T)                          298 - 848   
% % % S3         6 - 1017272.51          + 1231.20246     T       + 3621668.52     T^-1    848 - 1870  
% % % S3         6 - 110726059.     T^-2 - 207.930830     T ln(T)                          848 - 1870  
% % % S4         7 - 1015705.37          + 1229.35613     T       + 3621668.52     T^-1    848 - 1870  
% % % S4         7 - 110726059.     T^-2 - 207.930830     T ln(T)                          848 - 1870  
% % % S4         8 - 1022124.35          + 1274.89000     T       - 213.384000     T ln(T) 1870 - 2500 
% % % L1         9 - 1047181.60          + 1184.88110     T       + 3621668.52     T^-1    298 - 1870  
% % % L1         9 - 110726059.     T^-2 - 207.930830     T ln(T)                          298 - 1870  
% % % L1        10 - 1053600.58          + 1230.41497     T       - 213.384000     T ln(T) 1870 - 2500 
% % % ____________ _____________________ ________________________ ________________________ ___________ 
% % % 
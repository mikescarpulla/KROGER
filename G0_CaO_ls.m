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
G0_CaO_s = zeros(size(T));
G0_CaO_l = zeros(size(T));

% solid
masks1 = (T>298) .* (T<=2845);      %logical array same size as T with ones where true and 0 where false
masks2 = (T>2845) .* (T<=3500);
G0_CaO_s = G0_CaO_s + masks1.*(- 651262.658 + 573572.991*T.^(-1) - 17163131.3*T.^(-2) - 535.615998*T.^(0.5) + 376.676564*T - 58.7911706*T.*log(T));
G0_CaO_s = G0_CaO_s + masks2.*(- 676442.671 + 407.120846*T - 62.7600000*T.*log(T));

% liquid
maskl1 = (T>298) .* (T<=2845);      
maskl2 = (T>2845) .* (T<=3500);
G0_CaO_l = G0_CaO_l + maskl1.*(-571766.658 - 17163131.3*T.^(-2) + 573572.991*T.^(-1) + 348.735802*T - 535.615998*T.^(0.5) - 58.7911706*T.*log(T));
G0_CaO_l = G0_CaO_l + maskl2.*(-596946.671 + 379.180084*T - 62.7600000*T.*log(T));

% choose the min Go between s and l to choose the right phase at each T
G0_CaO_ls = min( cat(3, G0_CaO_s, G0_CaO_l),[],3);  % stack the two Go matrices along a dummy dimension, then take the min along that dimension.  Yes the [] is needed in the syntax for min() 



% now convert units to eV per Ga2O
G0_CaO_ls = G0_CaO_ls/(avo*q);   % eV/Ga2O molecule

% Now take Ptot and Xi into account.  
G0_CaO_ls = G0_CaO_ls + kB_eV*T.* log(X_i);

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen 
G0_CaO_ls(G0_CaO_ls==0) = Inf;

end


% % 
% % View Data  CaO     Units:  T(K) P(atm) Energy(J) Quantity(mol) 
% % Name: calcium oxide
% % 
% %   G(T) J/mol - 1 atm  
% % 
% %              G(T)                     G(T)                   G(T)                     T(K)        
% % ____________ ________________________ ______________________ ________________________ ___________ 
% % 
% % S1         1 - 651262.658             + 376.676564     T     + 573572.991     T^-1    298 - 2845  
% % S1         1 - 535.615998     T^0.5   - 17163131.3     T^-2  - 58.7911706     T ln(T) 298 - 2845  
% % S1         2 - 676442.671             + 407.120846     T     - 62.7600000     T ln(T) 2845 - 3500 
% % L1         3 - 571766.658             + 348.735802     T     + 573572.991     T^-1    298 - 2845  
% % L1         3 - 535.615998     T^0.5   - 17163131.3     T^-2  - 58.7911706     T ln(T) 298 - 2845  
% % L1         4 - 596946.671             + 379.180084     T     - 62.7600000     T ln(T) 2845 - 3500 
% % G1         5 78867.7875               + 205.751461     T     + 1.280263453E-02 T^ 2   298 - 1400  
% % G1         5 - 537947.001     T^-1    - 1.712449235E-06 T^ 3 - 10352.2747     ln(T)   298 - 1400  
% % G1         5 - 61.5491809     T ln(T)                                                 298 - 1400  
% % G1         6 37566279.6               - 406018.798     T     + 3.27569158     T^ 2    1400 - 3800 
% % G1         6 - 329377998.     T^-1    - 2.831966034E-05 T^ 3 - 14615548.9     ln(T)   1400 - 3800 
% % G1         6 - 1162.61296     T^ 1.5  + 5256160.29     T^0.5 + 48763.2209     T ln(T) 1400 - 3800 
% % G1         7 - 12975648.8             + 3153.92066     T     + 873451717.     T^-1    3800 - 6000 
% % G1         7 + 2367789.28     ln(T)   - 178046.562     T^0.5 - 282.261951     T ln(T) 3800 - 6000 
% % ____________ ________________________ ______________________ ________________________ ___________ 
% % 

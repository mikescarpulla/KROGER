function [G0_SnO_gv] = G0_SnO_gv(T, P_tot, X_i, P_units)
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
G0_SnO_gv = zeros(size(T));

% define masks needed based on Tranges given for the expressions for G0.  logical array same size as T with ones where true and 0 where false
mask1 = (T>298) .* (T<=1000); 
mask2 = (T>1000) .* (T<=2000); 

% solid
G0_SnO_gv = G0_SnO_gv + mask1.*(-298331.317 + 30.2261397*T + 140400*T.^(-1) - 221.320000*T.^(0.5) -37.96*T.*log(T));
G0_SnO_gv = G0_SnO_gv + mask2.*(9898.11022 + 2.11008105*T -6.515000000E-04*T.^(2) -34.7*T.*log(T));


% now convert units to eV per Ga2O
G0_SnO_gv = G0_SnO_gv/(avo*q);   % eV/Ga2O molecule

% Now take Ptot and Xi into account.  
% % % G0 = G0 + kB_eV*T.* log(X_i);  % condensed phase case (liq/solid)
G0_SnO_gv = G0_SnO_gv + kB_eV*T.* ( log(P_tot/P_ref) + log(X_i));  % gas/vapor case

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen 
G0_SnO_gv(G0_SnO_gv==0) = Inf;

end


% % % 
% % % View Data  SnO     Units:  T(K) P(atm) Energy(J) Quantity(mol) 
% % % Name: Tin Monoxide
% % % 
% % %   G(T) J/mol - 1 atm  
% % % 
% % %              G(T)                     G(T)                     G(T)                   T(K)        
% % % ____________ ________________________ ________________________ ______________________ ___________ 
% % % 
% % % S1         1              +      T       -  T^ 2 298 - 1273  
% % % S1         1    T ln(T)                                                 298 - 1273  
% % % G1         2 10508.1867               + 30.2261397     T       + 140400.000     T^-1  298 - 1000  
% % % G1         2 - 221.320000     T^0.5   - 37.9600000     T ln(T)                        298 - 1000  
% % % G1         3 9898.11022               + 2.11008105     T       - 6.515000000E-04 T^ 2 1000 - 2000 
% % % G1         3 - 34.7000000     T ln(T)                                                 1000 - 2000 
% % % ____________ ________________________ ________________________ ______________________ ___________ 
% % % 

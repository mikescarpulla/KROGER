function [G0_MgO_gv] = G0_MgO_gv(T, P_tot, X_i, P_units)
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
G0_MgO_gv = zeros(size(T));


% define masks needed based on Tranges given for the expressions for G0.  logical array same size as T with ones where true and 0 where false
mask1 = (T>298) .* (T<=1100);     
mask2 = (T>1100) .* (T<=4500); 
mask3 = (T>4500) .* (T<=6000); 

% gas/vapor
G0_MgO_gv = mask1.*(- 1102681.16 + 8518.40222*T + 10851502.7*T.^(-1) - 141244.711*T.^(0.5) + 5.878503055E-02*T.^(2) + 3.784973866E-06*T.^(3) + 435353.536*log(T) - 898.190474*T.*log(T));
G0_MgO_gv = G0_MgO_gv + mask2.*(- 2483196.14 + 70329447*T.^(-1) + 2388.43394*T + 3.419105015e-3*T.^(2) - 78828.9894*T.^(0.5) + 582033.945*log(T) - 243.395147*T.*log(T));
G0_MgO_gv = G0_MgO_gv + mask3.*(- 13298543.7 + 977136770*T.^(-1) + 3472.65363*T - 191013.773*T.^(0.5) + 2447324.28*log(T) - 305.763456*T.*log(T));



% now convert units to eV per Ga2O
G0_MgO_gv = G0_MgO_gv/(avo*q);   % eV/Ga2O molecule

% Now take Ptot and Xi into account.  
G0_MgO_gv = G0_MgO_gv + kB_eV*T.* ( log(P_tot/P_ref) + log(X_i));

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen 
G0_MgO_gv(G0_MgO_gv==0) = Inf;

end


%% for referecne, it's nice to copy the original data here in comments to allow proofreading.  For exaple this is the text file from FactSage for H2O.  
% % % 
% % % View Data  MgO     Units:  T(K) P(atm) Energy(J) Quantity(mol) 
% % % Name: magnesium oxide
% % % 
% % %   G(T) J/mol - 1 atm  
% % % 
% % %              G(T)                     G(T)                     G(T)                     T(K)        
% % % ____________ ________________________ ________________________ ________________________ ___________ 
% % % 
% % % S1         1 - 611541.379             + 420.064762     T       + 310577.002     T^-1    298 - 3098  
% % % S1         1 - 1184.79600     T^0.5   - 974102.005     T^-2    - 61.1096505     T ln(T) 298 - 3098  
% % % S1         2 - 662387.670             + 462.122764     T       - 66.9440000     T ln(T) 3098 - 3500 
% % % L1         3 - 554894.205             + 490.908663     T       + 1.571092000E-03 T^ 2   298 - 3098  
% % % L1         3 - 261375.798     T^-1    - 1184.79600     T^0.5   - 974102.005     T^-2    298 - 3098  
% % % L1         3 - 72.7955625     T ln(T)                                                   298 - 3098  
% % % L1         4 - 584985.686             + 437.137549     T       - 66.9440000     T ln(T) 3098 - 3500 
% % % G1         5 - 1102681.16             + 8518.40222     T       + 5.878503055E-02 T^ 2   298 - 1100  
% % % G1         5 + 10851502.7     T^-1    + 3.784973866E-06 T^ 3   + 435353.536     ln(T)   298 - 1100  
% % % G1         5 - 141244.711     T^0.5   - 898.190474     T ln(T)                          298 - 1100  
% % % G1         6 - 2483196.14             + 2388.43394     T       + 3.419105015E-03 T^ 2   1100 - 4500 
% % % G1         6 + 70329447.0     T^-1    + 582033.945     ln(T)   - 78828.9894     T^0.5   1100 - 4500 
% % % G1         6 - 243.395147     T ln(T)                                                   1100 - 4500 
% % % G1         7 - 13298543.7             + 3472.65363     T       + 977136770.     T^-1    4500 - 6000 
% % % G1         7 + 2447324.28     ln(T)   - 191013.773     T^0.5   - 305.763456     T ln(T) 4500 - 6000 
% % % ____________ ________________________ ________________________ ________________________ ___________ 
% % % 

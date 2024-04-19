function [G0_TiO_ls] = G0_TiO_ls(T, P_tot, X_i, P_units)
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
G0_TiO_s1 = zeros(size(T));
G0_TiO_s2 = zeros(size(T));
G0_TiO_liquid = zeros(size(T));
G0_TiO_ls = zeros(size(T));


% define masks needed based on Tranges given for the expressions for G0.  logical array same size as T with ones where true and 0 where false
mask1 = (T>298) .* (T<=2500);      % solid 1
mask2 = (T>298) .* (T<=2500) ;      % solid 2     
mask3 = (T>298) .* (T<=1200);      % liquid. 
mask4 = (T>1200) .* (T<=4500);

% solid1
G0_TiO_s1 = mask1.*(-558169.702 +255.475182*T -8.898664224E-03*T.^(2) +327010.589*T.^(-1) +1.105968351E-08*T.^(3) -41.9944927*T.*log(T));

% solid2
G0_TiO_s2 = mask2.*(-553971.140 +252.162502*T -8.899251973E-03*T.^(2) +327074.648*T.^(-1) +1.121961627E-08*T.^(3) -41.9953178*T.*log(T));

% liquid
G0_TiO_liquid =  mask3.*(-509468.783 +230.218069*T -8.858076521E-03*T.^(2) +326982.846*T.^(-1) -42.0139010*T.*log(T));
G0_TiO_liquid = G0_TiO_liquid + mask4.*(-526084.300 +410.418030*T -66.9440000*T.*log(T));

% choose the min Go, thus automatically selecting the right phase.  This is
% written assuming
G0_TiO_ls = min( cat(3,G0_TiO_s1, G0_TiO_s2, G0_TiO_liquid),[],3);  % stack the two Go matrices along a dummy dimension, then take the min along that dimension.  Yes the [] is needed in the syntax for min() 

% now convert units to eV per Ga2O
G0_TiO_ls = G0_TiO_ls/(avo*q);   % eV/Ga2O molecule

% Now take Ptot and Xi into account.  
G0_TiO_ls = G0_TiO_ls + kB_eV*T.* ( log(P_tot/P_ref) + log(X_i));

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen 
G0_TiO_ls(G0_TiO_ls==0) = Inf;

end


%% for referecne, it's nice to copy the original data here in comments to allow proofreading.  For exaple this is the text file from FactSage for H2O.  
% % % 
% % % View Data  TiO     Units:  T(K) P(atm) Energy(J) Quantity(mol) 
% % % Name: Titanium Oxide
% % % 
% % %   G(T) J/mol - 1 atm  
% % % 
% % %              G(T)                     G(T)                     G(T)                     T(K)        
% % % ____________ ________________________ ________________________ ________________________ ___________ 
% % % 
% % % S1         1 - 558169.702             + 255.475182     T       - 8.898664224E-03 T^ 2   298 - 2500  
% % % S1         1 + 327010.589     T^-1    + 1.105968351E-08 T^ 3   - 41.9944927     T ln(T) 298 - 2500  
% % % S1         2 - 565010.826             + 431.292509     T       - 66.9440000     T ln(T) 2500 - 2501 
% % % S2         3 - 553971.140             + 252.162502     T       - 8.899251973E-03 T^ 2   298 - 2500  
% % % S2         3 + 327074.648     T^-1    + 1.121961627E-08 T^ 3   - 41.9953178     T ln(T) 298 - 2500  
% % % S2         4 - 560811.474             + 427.972597     T       - 66.9440000     T ln(T) 2500 - 2501 
% % % L1         5 - 509468.783             + 230.218069     T       - 8.858076521E-03 T^ 2   298 - 1200  
% % % L1         5 + 326982.846     T^-1    - 42.0139010     T ln(T)                          298 - 1200  
% % % L1         6 - 526084.300             + 410.418030     T       - 66.9440000     T ln(T) 1200 - 4500 
% % % G1         7 78198.9955               + 125.670094     T       + 4.329207354E-03 T^ 2   298 - 1400  
% % % G1         7 - 423001.463     T^-1    - 3.836878386E-07 T^ 3   - 7615.16143     ln(T)   298 - 1400  
% % % G1         7 - 50.8736187     T ln(T)                                                   298 - 1400  
% % % G1         8 225627.296               - 728.339114     T       - 3.858717941E-03 T^ 2   1400 - 3700 
% % % G1         8 - 416191.214     T^-1    - 61603.3296     ln(T)   + 15770.0655     T^0.5   1400 - 3700 
% % % G1         8 + 34.3386587     T ln(T)                                                   1400 - 3700 
% % % G1         9 6685469.59               - 1273.18812     T       - 463897667.     T^-1    3700 - 6000 
% % % G1         9 - 1202306.56     ln(T)   + 86731.5508     T^0.5   + 56.8468382     T ln(T) 3700 - 6000 
% % % ____________ ________________________ ________________________ ________________________ ___________ 
% % % 

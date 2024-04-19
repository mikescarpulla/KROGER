function [G0_FeO_ls] = G0_FeO_ls(T, P_tot, X_i, P_units)
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
G0_FeO_s1 = zeros(size(T));
G0_FeO_l1 = zeros(size(T));
G0_FeO_ls = zeros(size(T));


% solid1
mask1 = (T>298) .* (T<=1644);      %logical array same size as T with ones where true and 0 where false
G0_FeO_s1 = G0_FeO_s1 + mask1.*(-322147.542  -330.687435*T  -1.530402997E-02*T.^(2)  +1266650.00*T.^(-1) +6003.60001*T.^(0.5) + 18.0244741*T.*log(T) );
mask2 = (T>1644) .* (T<=2000);
G0_FeO_s1 = G0_FeO_s1 + mask2.*(-299283.753 + 417.258468*T -68.1992000*T.*log(T)  );

% liquid1
mask3 = (T>298) .* (T<=1644);      %logical array same size as T with ones where true and 0 where false
G0_FeO_l1 = G0_FeO_l1 + mask3.*(-290958.454  -349.657168*T     -1.530402997E-02*T.^(2) + 1266650.00*T.^(-1)    + 6003.60001*T.^(0.5) + 18.0244741*T.*log(T)  );
mask4 = (T>1644) .* (T<=2000);
G0_FeO_l1 = G0_FeO_l1 + mask4.*(-268094.665  + 398.288735*T     -68.1992000*T.*log(T)   );


% choose the min Go, thus automatically selecting the right phase for all T's and P's.  
G0_FeO_ls = min( cat(3, G0_FeO_s1, G0_FeO_l1)  ,[],3);  % stack the Go matrices along a dummy dimension, then take the min along that dimension.  Yes the [] is needed in the syntax for min() 

% now convert units to eV per Ga2O
G0_FeO_ls = G0_FeO_ls/(avo*q);   % eV/Ga2O molecule

% Now take Ptot and Xi into account.  
G0_FeO_ls = G0_FeO_ls + kB_eV*T.* log(X_i);

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen 
G0_FeO_ls(G0_FeO_ls==0) = Inf;

end


% % % 
% % % View Data  FeO     Units:  T(K) P(atm) Energy(J) Quantity(mol) 
% % % Name: iron  monoxide
% % % 
% % %   G(T) J/mol -1 atm  
% % % 
% % %   G(T)          G(T)        G(T)     *T(K)        
% % % ____________ ________________________ ______________________ ________________________ ___________ 
% % % 
% % % S1         1 -322147.542  -330.687435*T     -1.530402997E-02*T.^(2)   298 -1644  
% % % S1         1 + 1266650.00*T.^(-1)    + 6003.60001*T.^(0.5) + 18.0244741*T.*log(T)  ) 298 -1644  
% % % S1         2 -299283.753  + 417.258468*T     -68.1992000*T.*log(T)  ) 1644 -2000 
% % % L1         3 -290958.454  -349.657168*T     -1.530402997E-02*T.^(2)   298 -1644  
% % % L1         3 + 1266650.00*T.^(-1)    + 6003.60001*T.^(0.5) + 18.0244741*T.*log(T)  ) 298 -1644  
% % % L1         4 -268094.665  + 398.288735*T     -68.1992000*T.*log(T)  ) 1644 -2000 
% % % G1         5 264803.572    + 73.9875439*T     + 2.417114297E-03*T.^(2)   298 -1600  
% % % G1         5 -253638.108*T.^(-1)    -2.194756174E-07 T^ 3 -5386.03812     ln(T)   298 -1600  
% % % G1         5 -45.0863911*T.*log(T)  )                298 -1600  
% % % G1         6 -835030.307  + 1185.92621*T     + 1.368558597E-03*T.^(2)   1600 -3900 
% % % G1         6 + 31173911.5*T.^(-1)    + 253414.297     ln(T) -37141.2387*T.^(0.5)   1600 -3900 
% % % G1         6 -140.037464*T.*log(T)  )                1600 -3900 
% % % G1         7 2605868.55    -473.856883*T     -151092816.*T.^(-1)    3900 -6000 
% % % G1         7 -436463.463     ln(T)   + 32771.7268*T.^(0.5) -2.41189124*T.*log(T)  ) 3900 -6000 
% % % ____________ ________________________ ______________________ ________________________ ___________ 
% % % 

function [G0_Fe2O3_ls] = G0_Fe2O3_ls(T, P_tot, X_i, P_units)
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
G0_Fe2O3_s1 = zeros(size(T));
G0_Fe2O3_s2 = zeros(size(T));
G0_Fe2O3_s3 = zeros(size(T));
G0_Fe2O3_ls = zeros(size(T));


% solid1
mask1 = (T>298) .* (T<=2500);      %logical array same size as T with ones where true and 0 where false
G0_Fe2O3_s1 = G0_Fe2O3_s1 + mask1.*(-861183.055 + 828.050052*T + 1453820.00*T.^(-1)  -137.008930*T.*log(T)  );

% solid 2
mask2 = (T>298) .* (T<=956); 
G0_Fe2O3_s2 = G0_Fe2O3_s2 + mask2.*(-809611.086 + 921.403602*T -2.740395108E-03*T.^(2) + 2788391.03*T.^(-1) + 6.890191675E-06*T.^(3)   -6.496518766E-09*T.^(4) -87604647.8*T.^(-2)  -146.858400*T.*log(T) );
mask3 = (T>956) .* (T<=1873); 
G0_Fe2O3_s2 = G0_Fe2O3_s2 + mask3.*(-802879.999 + 912.363938*T + 2788391.03*T.^(-1) -87604647.8*T.^(-2) -146.858400*T.*log(T)  );

% solid 3
mask4 = (T>298) .* (T<=1873); 
G0_Fe2O3_s3 = G0_Fe2O3_s3 + mask4.*(-801591.449  + 911.016249*T + 2788391.03*T.^(-1) -87604647.8*T.^(-2)    -146.858400*T.*log(T)  );

% choose the min Go, thus automatically selecting the right phase for all T's and P's.  
G0_Fe2O3_ls = min( cat(3, G0_Fe2O3_s1, G0_Fe2O3_s2, G0_Fe2O3_s3) ,[],3);  % stack the Go matrices along a dummy dimension, then take the min along that dimension.  Yes the [] is needed in the syntax for min() 

% now convert units to eV per Ga2O
G0_Fe2O3_ls = G0_Fe2O3_ls/(avo*q);   % eV/Ga2O molecule

% Now take Ptot and Xi into account.  
G0_Fe2O3_ls = G0_Fe2O3_ls + kB_eV*T.* log(X_i);

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen 
G0_Fe2O3_ls(G0_Fe2O3_ls==0) = Inf;

end


% % % 
% % % View Data  Fe2O3     Units:  T(K) P(atm) Energy(J) Quantity(mol) 
% % % Name: Hematite
% % % 
% % %   G(T) J/mol -1 atm  
% % % 
% % %   G(T)          G(T)          G(T)   *T(K)       
% % % ____________ ________________________ ________________________ ______________________ __________ 
% % % 
% % % S1         1 -861183.055  + 828.050052*T       + 1453820.00*T.^(-1)  298 -2500 
% % % S1         1 -137.008930*T.*log(T)  )                298 -2500 
% % % S1         1 + G(magnet)                  
% % % 
% % % S2         2 -809611.086  + 921.403602*T       -2.740395108E-03*T.^(2) 298 -956  
% % % S2         2 + 2788391.03*T.^(-1)    + 6.890191675E-06*T.^(3)   -6.496518766E-09*T.^(4) 298 -956  
% % % S2         2 -87604647.8*T.^(-2)    -146.858400*T.*log(T)  )             298 -956  

% % % S2         3 -802879.999  + 912.363938*T       + 2788391.03*T.^(-1)  956 -1873 
% % % S2         3 -87604647.8*T.^(-2)    -146.858400*T.*log(T)  )             956 -1873 
% % % 
% % % S3         4 -801591.449  + 911.016249*T       + 2788391.03*T.^(-1)  298 -1873 
% % % S3         4 -87604647.8*T.^(-2)    -146.858400*T.*log(T)  )             298 -1873 
% % % ____________ ________________________ ________________________ ______________________ __________ 
% % % 

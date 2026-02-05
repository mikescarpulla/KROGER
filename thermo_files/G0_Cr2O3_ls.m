function [G0_Cr2O3_ls] = G0_Cr2O3_ls(T, P_tot, X_i, P_units)
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
G0_Cr2O3_solid = zeros(size(T));
G0_Cr2O3_liquid = zeros(size(T));
G0_Cr2O3_ls = zeros(size(T));


% define masks needed based on Tranges given for the expressions for G0.  logical array same size as T with ones where true and 0 where false
mask1 = (T>=115) .* (T<=140);      % solid
mask2 = (T>140) .* (T<=298);
mask3 = (T>298) .* (T<=3000);
mask4 = (T>298) .* (T<=1800);      % liquid.  
mask5 = (T>1800) .* (T<=4500);

% solid
G0_Cr2O3_solid = mask1.*(-1140518.33 +7.570140550e-02*T.^(2) -7.047015000e-04*T.^(3) + 6.577516667e-07*T.^(4) -0.541157333*T.^(1.5) -8.96023508*T);
G0_Cr2O3_solid = G0_Cr2O3_solid + mask2.*(-1124581.89  +1274.14253*T  -4.118145350e-02*T.^(2)  -6534.39910*T.^(0.5)   -171.344603*T.*log(T) );
G0_Cr2O3_solid = G0_Cr2O3_solid + mask3.*(-1169514.24 +739.717283*T -3.900000000e-03*T.^(2) +1000000.00*T.^(-1)    -5.000000000e-08*T.^(3)   -121.440000*T.*log(T) );

% liquid
G0_Cr2O3_liquid = mask4.*(-1126252.44 +914.842281*T -3.420837616e-03*T.^(2) +1931683.05*T.^(-1)    +22280.8848*log(T)   -5069.97190*T.^(0.5) -141.819281*T.*log(T));
G0_Cr2O3_liquid = G0_Cr2O3_liquid + mask5.*(-1102615.27  +980.900759*T -156.690800*T.*log(T));

% choose the min Go, thus automatically selecting the right phase.  This is
% written assuming
G0_Cr2O3_ls = min( cat(3,G0_Cr2O3_solid, G0_Cr2O3_liquid),[],3);  % stack the two Go matrices along a dummy dimension, then take the min along that dimension.  Yes the [] is needed in the syntax for min() 


% now convert units to eV per Ga2O
G0_Cr2O3_ls = G0_Cr2O3_ls/(avo*q);   % eV/Ga2O molecule

% Now take Ptot and Xi into account.  
G0_Cr2O3_ls = G0_Cr2O3_ls + kB_eV*T.*  log(X_i);

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen 
G0_Cr2O3_ls(G0_Cr2O3_ls==0) = Inf;

end


%% for referecne, it's nice to copy the original data here in comments to allow proofreading.  For exaple this is the text file from FactSage for H2O.  
% % % 
% % % View Data  Cr2O3     Units:  T(K) P(atm) Energy(J) Quantity(mol) 
% % % Name: Chromium sesquioxide
% % % 
% % %   G(T) J/mol - 1 atm  
% % % 
% % %              G(T)                     G(T)                     G(T)                     T(K)        
% % % ____________ ________________________ ________________________ ________________________ ___________ 
% % % 
% % % S1         1 - 1140518.33             + 7.570140550E-02 T^ 2   - 7.047015000E-04 T^ 3   115 - 140   
% % % S1         1 + 6.577516667E-07 T^ 4   - 0.541157333     T^ 1.5 - 8.96023508     T       115 - 140   
% % % S1         1 + G(magnet)                                                                            
% % % S1         2 - 1124581.89             + 1274.14253     T       - 4.118145350E-02 T^ 2   140 - 298   
% % % S1         2 - 6534.39910     T^0.5   - 171.344603     T ln(T)                          140 - 298   
% % % S1         2 + G(magnet)                                                                            
% % % S1         3 - 1169514.24             + 739.717283     T       - 3.900000000E-03 T^ 2   298 - 3000  
% % % S1         3 + 1000000.00     T^-1    - 5.000000000E-08 T^ 3   - 121.440000     T ln(T) 298 - 3000  
% % % S1         3 + G(magnet)                                                                            
% % % L1         4 - 1126252.44             + 914.842281     T       - 3.420837616E-03 T^ 2   298 - 1800  
% % % L1         4 + 1931683.05     T^-1    + 22280.8848     ln(T)   - 5069.97190     T^0.5   298 - 1800  
% % % L1         4 - 141.819281     T ln(T)                                                   298 - 1800  
% % % L1         5 - 1102615.27             + 980.900759     T       - 156.690800     T ln(T) 1800 - 4500 
% % % ____________ ________________________ ________________________ ________________________ ___________ 
% % % 

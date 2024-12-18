function [G0_substance] = G0_substance(T, P_tot, X_i, P_units)
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
G0_substance_s1 = zeros(size(T));
G0_substance_s2 = zeros(size(T));
G0_substance = zeros(size(T));


% define masks needed based on Tranges given for the expressions for G0.  logical array same size as T with ones where true and 0 where false
mask1 = (T>0) .* (T<=800);      % solid 1
mask2 = (T>273) .* (T<=500) .* (P>=3.14);      % solid 2    We dont know where the transition from solid 1 to solid 2 is necessarily (could be f(T,P) not jsut f(T) )
mask3 = (T>500) .* (T<=3000);      % liquid.  Tmelt at stp is 500 K for this example

% solid1
G0_substance_s1 = mask1.*(-303622.945 + 198.232091*T - 36.2460*T.*log(T));

% solid2 - high pressure phase
G0_substance_s2 = mask2.*(A + B*T -C*T.*log(T) + func(T,P_tot)  );

% choose the min Go, thus automatically selecting the right phase.  This is
% written assuming
G0_substance = min( cat(3,G0_substance_s1 G0_substance_s2),[],3);  % stack the two Go matrices along a dummy dimension, then take the min along that dimension.  Yes the [] is needed in the syntax for min() 

% liquid
G0_substance = G0_substance + mask3.*(-256638.942 - 1924378.83*T.^(-1) - 1118.62474*T - 0.760342801*T.^(2) +5.3188744e-4*T.^(3) - 2.059132015e-7*T.^(4) + 203.118982*T.*log(T));

% now convert units to eV per Ga2O
G0_substance = G0_substance/(avo*q);   % eV/Ga2O molecule

% Now take Ptot and Xi into account.  For solids and liquids, it's only Xi
% that matters while gasses/vapors have the P/Pref term too.  Pick one of
% the lines below and delete the other one.  

G0_substance = G0_substance + kB_eV*T.* log(X_i);  %condensed phases

G0_substance = G0_substance + kB_eV*T.* ( log(P_tot/P_ref) + log(X_i));   %gas/vapor

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen 
G0_substance(G0_substance==0) = Inf;
G0_substance(isnan(G0_substance)) = Inf;   % set NaN to Inf too - probably this will arise at 0 K.  

end


%% for referecne, it's nice to copy the original data here in comments to allow proofreading.  For exaple this is the text file from FactSage for H2O.  
% % % 
% % % View Data  H2O     Units:  T(K) P(atm) Energy(J) Quantity(mol) 
% % % Name: Water
% % % 
% % %   G(T) J/mol - 1 atm  
% % % 
% % %              G(T)                     G(T)                   G(T)                     T(K)        
% % % ____________ ________________________ ______________________ ________________________ ___________ 
% % % 
% % % S1         1 - 303622.945             + 198.232091     T     - 36.2460000     T ln(T) 250 - 273   
% % % L1         2 - 256638.942             - 1118.62474     T     - 0.760349801     T^ 2   298 - 500   
% % % L1         2 - 1924378.83     T^-1    + 5.318874400E-04 T^ 3 - 2.059132015E-07 T^ 4   298 - 500   
% % % L1         2 + 203.118982     T ln(T)                                                 298 - 500   
% % % G1         3 - 255475.808             - 15.1731427     T     - 7.474858139E-03 T^ 2   298 - 1100  
% % % G1         3 + 13999.6597     T^-1    + 9.205931575E-08 T^ 3 + 1107.27182     ln(T)   298 - 1100  
% % % G1         3 - 25.7816397     T ln(T)                                                 298 - 1100  
% % % G1         4 152152.281               + 164.817017     T     - 8.054003823E-05 T^ 2   1100 - 4000 
% % % G1         4 - 12075583.2     T^-1    - 83128.2757     ln(T) + 5947.37003     T^0.5   1100 - 4000 
% % % G1         4 - 53.1457895     T ln(T)                                                 1100 - 4000 
% % % G1         5 - 4469439.52             + 1472.55188     T     + 298264202.     T^-1    4000 - 6000 
% % % G1         5 + 778290.990     ln(T)   - 64372.3401     T^0.5 - 155.190827     T ln(T) 4000 - 6000 
% % % ____________ ________________________ ______________________ ________________________ ___________ 
% % % 

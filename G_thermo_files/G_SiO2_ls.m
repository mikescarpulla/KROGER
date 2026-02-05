function [G0_SiO2_ls] = G0_SiO2_ls(T, P_tot, X_i, P_units)
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
% FAQ4: why keep the X_i for condensed phases?  You can thus use it to make an ideal solution of multiple phases.  You could make a function that adds a few phases and adds excess entropy and excess enthalpy to represent an alloy.... 
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
G0_SiO2_1 = zeros(size(T));  % since T and P are same size, we can use either to preallocate memory for this. 
G0_SiO2_2 = zeros(size(T));  % since T and P are same size, we can use either to preallocate memory for this. 
G0_SiO2_ls = zeros(size(T));  % since T and P are same size, we can use either to preallocate memory for this. 


% most stable solid
mask1 = (T>298) .* (T<=373);      %logical array same size as T with ones where true and 0 where false
mask2 = (T>373) .* (T<=848);   
mask3 = (T>848) .* (T<=850);   
G0_SiO2_1 = G0_SiO2_1 + mask1.*(-935388.521 +536.025342*T +1773342.00*T.^(-1) -961.103996*T.^(0.5) -81928061.6*T.^(-2) -80.0119918*T.*log(T));
G0_SiO2_1 = G0_SiO2_1 + mask2.*(-935486.558 +537.075695*T -4.220010851E-03*T.^(2) +1773342.00*T.^(-1)    +7.535450247E-06*T.^(3)   -5.045870528E-09*T.^(4) -961.103996*T.^(0.5)   -81928061.6*T.^(-2)    -80.0119918*T.*log(T) );
G0_SiO2_1 = G0_SiO2_1 + mask3.*( -876152.929 -104.031919*T -4.184000000E-02*T.*log(T) );

% 2nd solid
mask4 = (T>298) .* (T<=1996);      % logical array same size as T with ones where true and 0 where false
mask5 = (T>1996) .* (T<=3000);
G0_SiO2_2 = G0_SiO2_2 + mask4.*(-933315.349 +533.278543*T +1773342.00*T.^(-1) -961.103996*T.^(0.5)   -81928061.6*T.^(-2)    -80.0119918*T.*log(T)  );
G0_SiO2_2 = G0_SiO2_2 + mask5.*(-964566.440 +571.627480*T -85.7720000*T.*log(T) );





% choose the min Go, thus automatically selecting the right phase for each (T,P)  

% zero_index = find(G0_SiO2_2==0);
% G0_SiO2_2(zero_index) = Inf;
G0_SiO2_ls = min( cat(3,G0_SiO2_1, G0_SiO2_2),[],3);  % stack the two Go matrices along a dummy dimension, then take the min along that dimension.  Yes the [] is needed in the syntax for min() 

% now convert units to eV per Ga2O
G0_SiO2_ls = G0_SiO2_ls/(avo*q);   % eV/Ga2O molecule

% Now take Ptot and Xi into account.  
G0_SiO2_ls = G0_SiO2_ls + kB_eV*T.* log(X_i);

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen 
G0_SiO2_ls(G0_SiO2_ls==0) = Inf;

end


%% for referecne, it's nice to copy the original data here in comments to allow proofreading.  For exaple this is the text file from FactSage for H2O.  
% % % 
% % % View Data  SiO2     Units:  T(K) P(atm) Energy(J) Quantity(mol) 
% % % Name: silicon dioxide
% % % 
% % %   G(T) J/mol - 1 atm  
% % % 
% % %              G(T)                     G(T)                     G(T)                      T(K)        
% % % ____________ ________________________ ________________________ _________________________ ___________ 
% % % 
% % % S1         1 - 935388.521             + 536.025342     T       + 1773342.00     T^-1     298 - 373   
% % % S1         1 - 961.103996     T^0.5   - 81928061.6     T^-2    - 80.0119918     T ln(T)  298 - 373   
% % % S1         2 - 935486.558             + 537.075695     T       - 4.220010851E-03 T^ 2    373 - 848   
% % % S1         2 + 1773342.00     T^-1    + 7.535450247E-06 T^ 3   - 5.045870528E-09 T^ 4    373 - 848   
% % % S1         2 - 961.103996     T^0.5   - 81928061.6     T^-2    - 80.0119918     T ln(T)  373 - 848   
% % % S1         3 - 876152.929             - 104.031919     T       - 4.184000000E-02 T ln(T) 848 - 850   
% % % S2         4 - 933315.349             + 533.278543     T       + 1773342.00     T^-1     298 - 1996  
% % % S2         4 - 961.103996     T^0.5   - 81928061.6     T^-2    - 80.0119918     T ln(T)  298 - 1996  
% % % S2         5 - 964566.440             + 571.627480     T       - 85.7720000     T ln(T)  1996 - 3000 
% % % S3         6 - 945992.949             + 500.419346     T       - 9.103853469E-02 T^ 2    298 - 390   
% % % S3         6 + 2979047.54     T^-1    + 2.056362013E-04 T^ 3   - 1.741827512E-07 T^ 4    298 - 390   
% % % S3         6 - 159707687.     T^-2    - 75.3726680     T ln(T)                           298 - 390   
% % % S3         7 - 902938.323             - 57.0782497     T       - 4.184000000E-02 T ln(T) 390 - 392   
% % % S4         8 - 944111.184             + 480.752826     T       + 2979047.54     T^-1     298 - 1991  
% % % S4         8 - 159707687.     T^-2    - 75.3726680     T ln(T)                           298 - 1991  
% % % S4         9 - 961947.905             + 569.440181     T       - 85.7720000     T ln(T)  1991 - 3000 
% % % S5        10 - 926547.679             + 571.802887     T       - 1.010499806E-02 T^ 2    298 - 535   
% % % S5        10 + 1227679.99     T^-1    + 2.091769810E-05 T^ 3   - 1.623763504E-08 T^ 4    298 - 535   
% % % S5        10 - 1498.77200     T^0.5   - 46678699.1     T^-2    - 83.5135981     T ln(T)  298 - 535   
% % % S5        11 - 894636.817             - 74.9540946     T       - 4.184000000E-02 T ln(T) 535 - 537   
% % % S6        12 - 924997.138             + 566.999695     T       + 1227679.99     T^-1     298 - 1996  
% % % S6        12 - 1498.77200     T^0.5   - 46678699.1     T^-2    - 83.5135981     T ln(T)  298 - 1996  
% % % S6        13 - 961789.835             + 569.349438     T       - 85.7720000     T ln(T)  1996 - 3000 
% % % S7        14 - 905834.718             + 514.416682     T       + 1550000.46     T^-1     298 - 3000  
% % % S7        14 - 97633361.1     T^-2    - 6688.99846     ln(T)   - 78.0000108     T ln(T)  298 - 3000  
% % % S8        15 - 993685.735             + 353.061581     T       - 3.500999656E-03 T^ 2    298 - 3000  
% % % S8        15 + 6294493.28     T^-1    - 298799663.     T^-2    + 17010.0102     ln(T)    298 - 3000  
% % % S8        15 - 58.1199858     T ln(T)                                                    298 - 3000  
% % % L1        16 - 915415.778             + 562.199392     T       + 1227679.99     T^-1     298 - 1996  
% % % L1        16 - 1498.77200     T^0.5   - 46678699.1     T^-2    - 83.5135981     T ln(T)  298 - 1996  
% % % L1        17 - 952208.475             + 564.549135     T       - 85.7720000     T ln(T)  1996 - 3000 
% % % G1        18 - 308897.418             + 344.751947     T       - 2224.00000     T^0.5    298 - 1000  
% % % G1        18 - 76.0300000     T ln(T)                                                    298 - 1000  
% % % G1        19 - 317759.945             + 259.906267     T       - 1124.00000     T^0.5    1000 - 2400 
% % % G1        19 - 67.5000000     T ln(T)                                                    1000 - 2400 
% % % G1        20 - 331016.145             + 194.888775     T       - 7.985000000E-05 T^ 2    2400 - 6000 
% % % G1        20 - 61.3600000     T ln(T)                                                    2400 - 6000 
% % % ____________ ________________________ ________________________ _________________________ ___________ 
% % % % % % 
% % % 
% % % 
% % % These are foramtted right for matlab
% % % 
% % % 
% % % S3   6 -945992.949 +500.419346*T -9.103853469E-02*T.^(2)    298 -390   
% % % S3   6 +2979047.54*T.^(-1)    +2.056362013E-04*T.^(3)   -1.741827512E-07*T.^(4)    298 -390   
% % % S3   6 -159707687.*T.^(-2)    -75.3726680*T.*log(T)   298 -390   
% % % 
% % % S3   7 -902938.323 -57.0782497*T -4.184000000E-02*T.*log(T) 390 -392   
% % % 
% % % S4   8 -944111.184 +480.752826*T +2979047.54*T.^(-1)     298 -1991  
% % % S4   8 -159707687.*T.^(-2)    -75.3726680*T.*log(T)   298 -1991  
% % % 
% % % S4   9 -961947.905 +569.440181*T -85.7720000*T.*log(T)  1991 -3000 
% % % 
% % % S5  10 -926547.679 +571.802887*T -1.010499806E-02*T.^(2)    298 -535   
% % % S5  10 +1227679.99*T.^(-1)    +2.091769810E-05*T.^(3)   -1.623763504E-08*T.^(4)    298 -535   
% % % S5  10 -1498.77200*T.^(0.5)   -46678699.1*T.^(-2)    -83.5135981*T.*log(T)  298 -535   
% % % 
% % % S5  11 -894636.817 -74.9540946*T -4.184000000E-02*T.*log(T) 535 -537   
% % % 
% % % S6  12 -924997.138 +566.999695*T +1227679.99*T.^(-1)     298 -1996  
% % % S6  12 -1498.77200*T.^(0.5)   -46678699.1*T.^(-2)    -83.5135981*T.*log(T)  298 -1996  
% % % 
% % % S6  13 -961789.835 +569.349438*T -85.7720000*T.*log(T)  1996 -3000 
% % % 
% % % S7  14 -905834.718 +514.416682*T +1550000.46*T.^(-1)     298 -3000  
% % % S7  14 -97633361.1*T.^(-2)    -6688.99846*log(T)   -78.0000108*T.*log(T)  298 -3000  
% % % 
% % % S8  15 -993685.735 +353.061581*T -3.500999656E-03*T.^(2)    298 -3000  
% % % S8  15 +6294493.28*T.^(-1)    -298799663.*T.^(-2)    +17010.0102*log(T)    298 -3000  
% % % S8  15 -58.1199858*T.*log(T)    298 -3000  
% % % 
% % % L1  16 -915415.778 +562.199392*T +1227679.99*T.^(-1)     298 -1996  
% % % L1  16 -1498.77200*T.^(0.5)   -46678699.1*T.^(-2)    -83.5135981*T.*log(T)  298 -1996  
% % % 
% % % L1  17 -952208.475 +564.549135*T -85.7720000*T.*log(T)  1996 -3000 
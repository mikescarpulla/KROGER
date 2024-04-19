function [G0_Ga_gv] = G0_Ga_gv(T, P_tot, X_i, P_units)
% gives the thermodynamic variables of Ga vapor as function of T, P based
% on Matvei Zinkevich and Fritz Aldinger J. Am. Ceram. Soc., 87 [4] 683–91
% (2004)  
% T in K, P in atm or Torr (must put in string argument 'atm' or 'Torr' etc
%
% G computed in kJ/mol then after these are computed we
% convert to eV
% T and P should be vectors

% constants
q = 1.602176634e-19;
avo = 6.0221409e+23;
kB_eV = 8.617333262e-5;


if strcmp(P_units,'atm')
    P_ref = 1;
elseif strcmp(P_units,'Torr')
    P_ref = 760;   
elseif strcmp(P_units,'Bar')
    P_ref = 1;
elseif strcmp(P_units,'Pa')
    P_ref = 1e5;
else
    error('Units of pressure must be atm, Torr, Pa, or Bar unless the ln(P/Pref) term is changed')
end


T = T(:);  % make T col vector
nT=size(T,1);
P_tot = P_tot(:)';  % make P a row vector
nP = size(P_tot,1);

% note here T and P are vectors
T = T*ones(1,nP);   % copy T to make a matrix with as many columns as P has entries
P_tot = ones(nT,1)*P_tot;   % similarly for P - copy down as many rows as T has entries
G0_Ga_gv = zeros(size(T));

mask1 = (T>0) .* (T<=600);  % 0 to 600     logical array same size as T with ones where true and 0 where false
mask2 = (T>600) .* (T<=1400);   % 600-1400
mask3 = (T>1400) .* (T<=6000);   % 1400 - 6000


% below Tmelt (zeros except for below Tmelt entires)
% T region 1
G0_Ga_gv = G0_Ga_gv + mask1.*( 259072.278 + 88.0130706*T - 38.71057*T.*log(T) + 0.01053784*(T.^2) - 9.86907833E-7*(T.^3) + 338489.2*(T.^-1) );

% T region 2
G0_Ga_gv = G0_Ga_gv + mask2.*(263812.519 + 33.4871435*T - 30.75007*T.*log(T) + 0.00537745*(T.^2) - 5.46534E-7*(T.^3) - 150942.65*(T.^-1) );

% T region 3
G0_Ga_gv = G0_Ga_gv + mask3.*(270292.501 - 28.1810494*T - 21.9834*T.*log(T) + 3.192416E-4*(T.^2) - 1.46299133E-8*(T.^3) - 992093*(T.^-1) );


% now convert units to eV per Ga atom
G0_Ga_gv = G0_Ga_gv/(avo*q);   % eV/Ga atom


% Now add in the pressure (zero for condensed phases) 
G0_Ga_gv = G0_Ga_gv + kB_eV*T.*( log(P_tot/P_ref) + log(X_i));

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen 
G0_Ga_gv(G0_Ga_gv==0) = Inf;


end

% % % View Data  Ga     Units:  T(K) P(atm) Energy(J) Quantity(mol) 
% % % Name: Gallium
% % % 
% % %   G(T) J/mol - 1 atm  
% % % 
% % %              G(T)                     G(T)                   G(T)                     T(K)        
% % % ____________ ________________________ ______________________ ________________________ ___________ 
% % % 
% % % S1         1 - 78802.3587             - 158.582219     T     - 1.902578477E-02 T^ 2   298 - 500   
% % % S1         1 + 1096609.41     T^-1    + 15690.5714     ln(T) + 13.2583508     T ln(T) 298 - 500   
% % % L1         2 - 2295.65181             + 118.648473     T     + 7.412438567E-06 T^ 2   298 - 3000  
% % % L1         2 + 68793.1928     T^-1    - 15529672.3     T^-2  - 26.6221078     T ln(T) 298 - 3000  
% % % G1         3 288077.538               + 237.110940     T     + 2.425536013E-02 T^ 2   298 - 1000  
% % % G1         3 - 52295.0586     T^-1    - 3.215573323E-06 T^ 3 - 6702.90489     ln(T)   298 - 1000  
% % % G1         3 - 59.3878806     T ln(T)                                                 298 - 1000  
% % % G1         4 251052.221               - 33.4923382     T     - 4.572503877E-05 T^ 2   1000 - 2477 
% % % G1         4 - 590744.590     T^-1    + 4193.45617     ln(T) - 361.608560     T^0.5   1000 - 2477 
% % % G1         4 - 20.9646567     T ln(T)                                                 1000 - 2477 
% % % G1         5 922873.604               - 645.290262     T     - 1.083065439E-03 T^ 2   2477 - 4000 
% % % G1         5 - 24222767.7     T^-1    - 146967.288     ln(T) + 19865.7011     T^0.5   2477 - 4000 
% % % G1         5 + 32.4569489     T ln(T)                                                 2477 - 4000 
% % % G1         6 - 3636483.14             + 1134.15387     T     + 248935413.     T^-1    4000 - 6000 
% % % G1         6 + 736473.569     ln(T)   - 61819.1627     T^0.5 - 112.379908     T ln(T) 4000 - 6000 
% % % ____________ ________________________ ______________________ ________________________ ___________ 
% % % 






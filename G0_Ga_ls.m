function [G0_Ga_ls] = G0_Ga_ls(T, P_tot, X_i, P_units)
% gives the thermodynamic variables of Ga solid nad liquid as function of T, P based
% on Matvei Zinkevich and Fritz Aldinger J. Am. Ceram. Soc., 87 [4] 683–91
% (2004)  (where the ZA comes from)
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


T = T*ones(1,nP);   % copy T to make a matrix with as many columns as P has entries
P_tot = ones(nT,1)*P_tot;   % similarly for P - copy down as many rows as T has entries
G0_Ga_s = zeros(size(T));
G0_Ga_l = zeros(size(T));

mask1 = (T>0) .* (T<=302.92);  % 0 to Tm     logical array same size as T with ones where true and 0 where false
mask2 = (T>302.92) .* (T<=4000);   % above Tm


% below Tmelt (zeros except for below Tmelt entires)
% solid phase is more stable below Tm
G0_Ga_s = G0_Ga_s + mask1.*( -21312.331 + 585.263691*T -108.228783*T.*log(T) + 0.227155636*T.^(2) -1.18575257E-4*T.^(3) + 439954*T.^(-1) ); 
G0_Ga_l = G0_Ga_s + mask1.*( 5491.31  -18.073718*T  -7.0154E-17*T.^(7));

G0_Ga_s = G0_Ga_s + mask2.*( -7055.646 +132.7302*T -26.0692906*T.*log(T) +1.506E-4*T.^(2) -4.0173E-8*T.^(3) -118332*T.^(-1) + 1.64554E-23.*T.^(-9)  );
G0_Ga_l = G0_Ga_s + mask2.*( 5666.446 -18.680788*T  -1.64554E-23*T.^(-9)  );

G0_Ga_ls = min([G0_Ga_l G0_Ga_s],[],2);  %Take the min value as the stable phase vs T.  The way these are written, we need the solid to calc the liquid even when the iq is the stable phase


% now convert units to eV per Ga atom
G0_Ga_ls = G0_Ga_ls/(avo*q);   % eV/Ga atom

% Now take Ptot and Xi into account.  
G0_Ga_ls = G0_Ga_ls + kB_eV*T.* log(X_i);

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen 
G0_Ga_ls(G0_Ga_ls==0) = Inf;

end






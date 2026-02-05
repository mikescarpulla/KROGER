function [G0_GaO_gv] = G0_GaO_gv(T,P_tot, X_i, P_units)
% gives the thermodynamic variables of GaO as function of T, P based
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
P_tot = P_tot*ones(nT,1);   % similarly for P - copy down as many rows as T has entries
G0_GaO_gv = zeros(size(T));

mask1 = (T>0) .* (T<=800);  % 0 to 800K     logical array same size as T with ones where true and 0 where false
mask2 = (T>800) .* (T<=1500);   % 800-1500 K
mask3 = (T>1500) .* (T<=4000);   % up to 4000 K


% up to 800 K 
G0_GaO_gv = G0_GaO_gv + mask1.*(136904.191 -22.25862*T -30.49045*T.*log(T) -0.0048738965*(T.^2) - 2.51268E-7*(T.^3) + 57767.3 .*(T.^-1)  );

% up to 1500 K
G0_GaO_gv = G0_GaO_gv + mask2.*(137661.642 -54.439509*T -25.24208*T.*log(T) -0.01210693*(T.^2) + 1.273842E-6*(T.^3) + 207293.6 .*(T.^-1) );

% up to 4000  K
G0_GaO_gv = G0_GaO_gv + mask3.*(109485.47 +175.498025*T -57.18317*T.*log(T) +0.0036644975*(T.^2) -1.63582983e-7*(T.^3) + 4743017 .*(T.^-1) );


% now convert units to eV per O2
G0_GaO_gv = G0_GaO_gv/(avo*q);   % eV/O2 molecule

% Now add in the pressure  GaO is a gas under 
G0_GaO_gv = G0_GaO_gv + kB_eV*T.*(log(P_tot/P_ref) + log(X_i));

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen 
G0_GaO_gv(G0_GaO_gv==0) = Inf;

end






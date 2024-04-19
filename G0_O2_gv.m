function [G0_O2_gv] = G0_O2_gv(T, P_tot, X_i, P_units)
% gives the thermodynamic variables of O from O2 as function of T, P based
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
G0_O2_gv = zeros(size(T));

mask1 = (T>0) .* (T<=900);  % 0 to 900K     logical array same size as T with ones where true and 0 where false
mask2 = (T>900) .* (T<=3700);   % 900-3700 K 

% up to 900 K 
G0_O2_gv = G0_O2_gv + mask1.*(-6960.6927  -51.1831467*T -22.25862*T.*log(T) -0.01023867*(T.^2) + 1.339947e-6*(T.^3) -76749.55.*(T.^-1)  );

% up to 3700 K
G0_O2_gv = G0_O2_gv + mask2.*(-13136.0174 + 24.7432966*T -33.55726*T.*log(T) -0.0012348985*(T.^2) + 1.66943333e-8*(T.^3) - 539886.*(T.^-1) );


% now convert units to eV per O2
G0_O2_gv = G0_O2_gv/(avo*q);   % eV/O2 molecule

% Now add in the pressure  
G0_O2_gv = G0_O2_gv + kB_eV*T.*( log(P_tot/P_ref) + log(X_i));

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen 
G0_O2_gv(G0_O2_gv==0) = Inf;


end






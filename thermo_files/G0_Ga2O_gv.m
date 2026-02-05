function [G0_Ga2O_gv] = G0_Ga2O_gv(T,P_Ga2O,P_units)
% gives the thermodynamic variables of Ga2O as function of T, P based
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
P_Ga2O = ones(nT,1)*P_Ga2O;   % similarly for P - copy down as many rows as T has entries
G0_Ga2O_gv = zeros(size(T));

mask1 = (T>0) .* (T<=900);  % 0 to 900K     logical array same size as T with ones where true and 0 where false
mask2 = (T>900) .* (T<=6000);   % 900-6000 K 
% up to 900 K 
G0_Ga2O_gv = G0_Ga2O_gv + mask1.*(-116546.827  +58.7325575*T -50.24181*T.*log(T) -0.006764885*(T.^2) + 1.088093e-6*(T.^3) -233385.05.*(T.^-1)  );
% up to 6000 K
G0_Ga2O_gv = G0_Ga2O_gv + mask2.*(-120662.748 + 111.155971*T -58.11384*T.*log(T) -1.699258e-5*(T.^2) + 5.77646167e-10*(T.^3) + 643465.*(T.^-1) );


% now convert units to eV per Ga2O
G0_Ga2O_gv = G0_Ga2O_gv/(avo*q);   % eV/Ga2O molecule

% Now add in the pressure and mole fraction  
G0_Ga2O_gv = G0_Ga2O_gv + kB_eV*T.*(log(P_Ga2O/P_ref) + log(X_Ga2O));

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen 
G0_Ga2O_gv(G0_Ga2O_gv==0) = Inf;

end
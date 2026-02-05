function [G0_Ga2O3_ls] = G0_Ga2O3_ls(T, P_tot, X_i, P_units)
% gives the Gibbs energy of beta-Ga2O3 as function of T and P
% T vector in K
% G comes out in eV/Formula unit
% values from Zinkevich and Aldinger JACS (2004)
% convert to eV/FU at the end
% the only thing P is used for here is making the output the right size


% now convert units to eV per FU of b-Ga2O3 
q = 1.602176634e-19;  % J/eV
avo = 6.0221409e+23;  % formula units/mol
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
G0_Ga2O3_ls = zeros(size(T));  %preallocate space


G0_Ga2O3_ls = -1127917 + 684.8332*T - 112.3935*T.*log(T) - 0.00796268819*(T.^2) + 1080114./T ;   % in J/mol


G0_Ga2O3_ls = G0_Ga2O3_ls/(q*avo);    % this should be in eV/FU


% Now take Ptot and Xi into account.  
G0_Ga2O3_ls = G0_Ga2O3_ls + kB_eV*T.*log(X_i);

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen 
G0_Ga2O3_ls(G0_Ga2O3_ls==0) = Inf;


end






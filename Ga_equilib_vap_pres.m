function [Ga_equilib_vap_pres] = Ga_equilib_vap_pres(T, P_units) 
%% Gives the equilibrium vapor pressure of Ga over Ga vs T 
% T is a vector

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
    error('Units of pressure must be atm, Torr, Pa, or Bar')
end

T = T(:);   %make T a col vactor

% compute the pressure-independent difference in G0 between the vapor and
% condensed phases.  
[G0_Ga_ls] = G0_Ga_ls(T, P_ref, 1, P_units);   
[G0_Ga_gv] = G0_Ga_gv(T, P_ref, 1, P_units);

Ga_equilib_vap_pres = exp(-(G0_Ga_gv-G0_Ga_ls)./(kB_eV*T));

end


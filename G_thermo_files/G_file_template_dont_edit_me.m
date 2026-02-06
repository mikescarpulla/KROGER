function [G_substance] = G_substance(T, P_tot, X_i, P_units)
%
%  each substance should have a function G_substance_ls.m for condensed
%  phases, and G_substance_gv.m for gas/vapor.  See FAQs below.  
% 
% this function gives the standard Gibbs energy Go of substance which is a
% function of T only.  Then it turns it into a function of T, Ptot, and X  (mole fraction).  
% for liquids/solids we only add in ln(X) term, while for gasses we add
% ln(X) + ln(Ptot/pref) and since x = Ppartial/Ptot this ends up being the
% same as ln(Ppartial/Pref)
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

G_substance = G0_substance + kB_eV*T.* log(X_i);  %condensed phases

G_substance = G0_substance + kB_eV*T.* ( log(P_tot/P_ref) + log(X_i));   %gas/vapor

% set any that are zero becasue of masking to infintiy so it produces an
% obvious error that can be seen 
G_substance(G_substance==0) = Inf;
G_substance(isnan(G_substance)) = Inf;   % set NaN to Inf too - probably this will arise at 0 K.  

end

% Inverts the Fermi-Dirac integral to find eta, assuming a 3D density of states dos = c*sqrt(E-Eg).  Uses root-finding on
% residual function for forward integral.  eta = (Ef-Ec)/kBT.
% this means eta is positive if Ef is above Ec and negative if it is in the
% gap
% n, Nc3D must be scalars.  Note this will work for holes too - up
% to user to change eta into Ef relatie to the right band.


function [eta] = eta_Fermi_Dirac(n, Nc3D)

if ~isscalar(n) || ~isscalar(Nc3D)
    error('n and Nc3D must be scalars')
elseif n<0 || ~isreal(n)
    error('n must be a non-negative real number');
elseif Nc3D<0 || ~isreal(Nc3D)
    error('Nc3D must be a non-negative real number');
end

% sz_n = size(n);
% if ~all(sz_n==size(Nc3D))
%     error('n and Nc3D have to be same shape')
% elseif min(sz_n)~=1
%     error('n and Nc3D have to be vectors')
% end
%
% % Initialize output
% eta = zeros(sz_n);

% Set options for fzero outside the for loop
options = optimset('TolX', 1e-12, 'Display', 'off');


% for i = 1:max(sz_n)



% Highly Degenerate - n>>NCD
if n <= Nc3D
    eta_guess = log(n / Nc3D);
elseif n > Nc3D
    eta_guess = (3/4*sqrt(pi) * n/Nc3D)^(2/3);
end

% Define residual function
residual = @(eta_test) n_Fermi_Dirac(eta_test, Nc3D) - n;

% Find the root
try
    eta = fzero(residual, eta_guess, options);
catch
    % If single guess fails
    % eta_lower = min(eta_guess - 10, -20);
    % eta_upper = max(eta_guess + 10, 50);
    eta_lower = eta_guess - 20;
    eta_upper = eta_guess + 20;
    try
        eta = fzero(residual, [eta_lower, eta_upper], options);
    catch
        eta_lower = eta_lower - 20;
        eta_upper = eta_upper + 20;
        eta = fzero(residual, [eta_lower, eta_upper], options);
    end
end

% end % for loop

end  %function

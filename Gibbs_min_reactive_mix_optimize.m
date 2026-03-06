function [out] = Gibbs_min_reactive_mix_optimize(in)


if strcmp(in.G_units,"eV/molecule")
    in.RkB = 8.617333262e-5;
elseif strcmp(in.G_units,"J/mol")
    in.RkB = 8.31446261815324;
else
    error('Must choose eV/molecule or J/mole as units of G')
end

if strcmp(in.P_units,'atm')
    in.P_ref = 1;
elseif strcmp(in.P_units,'Torr')
    in.P_ref = 760;
elseif strcmp(in.P_units,'Bar')
    in.P_ref = 1;
elseif strcmp(in.P_units,'Pa')
    in.P_ref = 1e5;
else
    error('Units of pressure must be atm, Torr, Pa, or Bar unless the log(P/Pref) term is changed')
end


% Default options for max iterations and convergance tolerance ----
% if ~isfield(in,'max_iter'); in.max_iter = 50; end
% if ~isfield(in,'converge_tol'); in.converge_tol = 5e-5; end
% if ~isfield(in,'min_size'); in.min_size = 1e-6; end

% copy over all the problem setup info
out = in;

% # elements and # species from a_ij matrix
[NE, NS] = size(in.a_ij);

% rotate the species moles to a column vector so it can multiply the a_ij
% without transposing it, which is needed for fmincon
in.init_species_moles = in.init_species_moles(:);

% set any moles that are zero to a number much smaller than the smallest
% nonzero one
zero_moles_index = in.init_species_moles==0;
nonzero_moles_index = ~zero_moles_index;
in.init_species_moles(zero_moles_index) = min(in.init_species_moles(nonzero_moles_index),[],"all")/in.max_min_mole_ratio;

% Compute initial total moles of elements.  This is conserved throguh calc.
out.init_element_moles = in.a_ij * in.init_species_moles;  % (NE x 1).  This is the number of moles of each element i.  It will be always conserved.    The (:) makes sure the flows/moles are used as a column vector
% out.init_gas_moles = in.species_states * in.init_species_moles; 
out.init_G_tot = Total_G_const_P(in.init_species_moles);

% Define inequalities
A = [];
b = [];

% Define linear equality constraints as Ax=b where x is the species_moles col vector
% element conservation
Aeq1 = in.a_ij;
beq1 = out.init_element_moles; %col vec with numel=NE

% % num gas moles conservation - this was a mistake since the flow velocity
% can change in our problem.  If it was a closed reactor with constant
% pressure, then gas moles would be constrained but that means only
% reactions that conserve total gas moles would be possible
% Aeq2 = in.species_states;  % here we add a row to the a_ij matrix that has a 1 for each species if it is a gas and a zero if it's condensed.
% beq2 = out.init_gas_moles; %scalar

% build up total matrix and vector of linear constraints 
% Aeq = [Aeq1; Aeq2];
% beq = [beq1; beq2];
Aeq = Aeq1;
beq = beq1;

% lb = in.min_size * ones(1,NS);  % makes sure numbes of moles can't go to zero or below some minimum
lb = [0 0 0];  % if this doesnt work, set some min number of moles as a very small fraction of the initial total moles (to get scaling right - as opposed to a hard coded min number since the total num moles could be small)
ub = [];
nonlcon = [];

mole_step = min(in.init_species_moles)/in.converge_factor;
G_step = abs(out.init_G_tot/in.converge_factor);

% options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
% options = optimoptions('fmincon','Display','iter','MaxIterations',in.max_iter);
options = optimoptions('fmincon','Display','iter', 'StepTolerance', mole_step, 'OptimalityTolerance', G_step);

% Determine equilibrium mole numbers and equilibrium G for the mixture
[out.equilib_species_moles, out.equilib_G_tot] = fmincon(@(species_moles) Total_G_const_P(species_moles), in.init_species_moles, A, b, Aeq, beq, lb, ub, nonlcon, options);

out.equilib_element_moles = in.a_ij * out.equilib_species_moles;
out.rel_element_mole_err = (out.equilib_element_moles - out.init_element_moles)./out.init_element_moles;




%% subroutines
    function G_tot = Total_G_const_P(species_moles)
        % % Total Gibbs Energy of Current Mixture.  Will come out in units defined.
        %  As cof now, the condensed vs gas phases have different G0 file names,
        %  with _gv taking pressure into account

        tot_moles = sum(species_moles);  % this sums over all species regardless of state.  In our G0 functions we multiply by mole fraction X for all species in a mix, and by pressure only for the gases
        species_X = species_moles/tot_moles; % this is the mole fraction of the total mixture, regardless of state
        g_vec = zeros(1,NS);  %rowvec
        for w = 1:NS   % loop over each species, compute the gibbs energy for that species.  Units dont matter (J/mol, eV/formula unit, etc).  We just want to minimize it.
            g_vec(w) = feval(in.Gibbs_func_names(w), in.T, in.P_tot, species_X(w), in.P_units); % this calls the G0 function for this species, accounting for total pressure and mole fractions
        end
        G_tot = g_vec * species_moles;
    end





end  %entire function
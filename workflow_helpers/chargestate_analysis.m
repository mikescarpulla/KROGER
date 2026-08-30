function [CTLVec,EVec,qPair,CTL,CTLMinVec,EMinVec,qMinPair,minCTL,qMinVec,E0MinVec] = chargestate_analysis(qVec,E0Vec,T,numKbt)
% chargestate_analysis calculates and analyzes the charge transition levels 
% for a defect with various charge states. The code finds intersections 
% between charge states, and identifies minimum energy points for each Fermi energy 
% level, then organizes the results for easy access. This function also
% produces 3 plots: 
%   1.) All charge state lines with intersections labeled ([q1,q2],CTL)
%   2.) Minimum charge state lines with intersection labeled ([q1,q2],CTL)
%   3.) Lowest charge state piecewise function with intersections labeled
%   ([q1,q2],CTL)
% Inputs: ***qVec and E0Vec MUST have matching indices***
%   - qVec: a column vector containing charge state values for a defect
%   - E0Vec: a column vector containing baseline energy values for a defect
%   (MUST BE IN eV)
%   - T: the temperature (MUST BE IN KELVIN)
%   - numKbt: how many kbt to consider (for accuracy) [#]
% Outputs: 
%   - CTLVec: a column vector of ALL fermi energies at charge transition levels (intersections).
%   - EVec: a column vector of ALL formation energies at charge transition
%   levels (intersections). 
%   - qPair: an N x 2 array that specifies the indexes (respectively) of the charge states
%   ([q1,q2])
%   - CTL: a N x 4 array that specifies respectively (by column):
%   q1, q2, Ef, E (a summary). 
%   - CTLMinVec: a column vector of fermi energies for MINIMUM formation energies of charge
%   transition levels (intersections). 
%   - EMinVec: a column vector of MINIMUM formation energies for charge
%   transition level's (intersections). 
%   - qMinPair: an Nmin x 2 array that specifies the indexes (respectively)
%   of the charge states ([q1,q2]). 
%   - minCTL: a Nmin x 4 array that specifies respectively (by column):
%   q1min, q2min, Efmin, Emin (a summary). 
%   - qMinVec: a column vector containing the unique values of q which gives
%   the lowest intersection points (indexes match E0MinVec).
%   - E0MinVec: a column vector containing the unique values of E0 which
%   gives the lowest intersection points (indexes match qMinVec).

%% Preperation of Variables and Errors
    % Covert column vectors to row vectors
    qVec = qVec'; % a [#]
    E0Vec = E0Vec'; % in [eV]

    % Define Variables
    % kb = 8.617333262145e-5; % Stephan Boltzmanns Constant in [eV/K]
    numq = numel(qVec); % Number of charge states [#]
    numE0 = numel(E0Vec); % Number of baseline energies [#]
    CTLVec = []; % an empty vector to store intersection CTL (x) values
    CTLIdx = 1; % an index counter to populate CTLvec 
    EVec = []; % an empty vector to store intersection E (y) values
    EIdx = 1; % an index counter to populate EVec 
    qPair = []; % Matrix to store pairs of charge states causing each intersection
    tolerance = 1e-6; % Set a small tolerance that depends on temperature (T) and numkbt


    % Throwing up Errors if vectors don't match
    if size(qVec, 1) ~= 1 || size(E0Vec, 1) ~= 1
        error('Charge state (q) and energy (E0) vectors must be row vectors.');
    elseif numq ~= numE0
        error('Charge state vector and baseline energy vector must have the same length.');
    end

%% Calculate Charge Transition Levels
    % Loop through all unique pairs of charge states
    for i = 1:numq-1
        for j = i+1:numq
        
        % Solve Ef and stored in CTLVec(CTLIndex)
        CTLVec(CTLIdx) = -(E0Vec(j)-E0Vec(i))/((qVec(j)-qVec(i))); 
        
        % Plug CTLVec back in to find E (formation energy)
        EVec(CTLIdx) = E0Vec(i) + qVec(i) * CTLVec(CTLIdx);

        % update CTL index and EIdx for next loop
        CTLIdx = CTLIdx + 1;

        % Store the charge state pair that caused this intersection
        qPair(:, end+1) = [qVec(i); qVec(j)]; 
        % ?????????????

       end
    end

%% Check Lowest Formation Energy 
    % Make new vectors which are copys of the original vectors in order to
    % specify only the lowest formation CTLs
    EMinVec = EVec; % copy of EVec
    CTLMinVec = CTLVec; % copy of CTLVec
    qMinPair = qPair; % copy of qMinPair

    % Check to see that the charge transition level is the lowest energy
    %   point. Reconfiguring in EMinVec and EfMinVec.
    for p = length(CTLMinVec):-1:1 % Iterate from end to beginning
        % Define Current Ef and E values
        EfCurr = CTLVec(p); % Current Ef value
        EComp = EVec(p); % E value to compare
        % Loop through all the charge states
        for z=1:length(qVec)
            % Calculate E from the q and E0 Value
            ECurr = E0Vec(z) + qVec(z)*EfCurr; 
            % If statement to replace EComp if ECurr is smaller
            if ECurr < EComp - tolerance
                CTLMinVec(p) = []; % deleting row
                EMinVec(p) = []; % deleting row
                qMinPair(:,p) = []; % deleting row
                break; % if we find a lower value
            end
        end
    end

%% Transpose CTLVec, EVec, qPair, EMinVec, CTLMinVec, and qMinPair
    % Transpose output vectors to column vectors
    CTLVec = CTLVec(:); 
    EVec = EVec(:); 
    qPair = qPair';
    EMinVec = EMinVec(:); 
    CTLMinVec = CTLMinVec(:);
    qMinPair = qMinPair';

%% Concatenate and Sort CTL and CTLMin Vectors 
    % Sort CTL by the Ef values in the third column
    CTL = sortrows([qPair,CTLVec,EVec],3);
    % Sort minCTL by the Ef values in the third column
    minCTL = sortrows([qMinPair,CTLMinVec,EMinVec],3);

%% Calculating qMinVec and E0MinVec
    % This section takes in an nx2 array of the minimum intersection charge state
    %   (q) values, and spits out a unique row vector of each of those q values and
    %   the corresponding E0 Values.
    % Flatten the input array into a column vector
    flatArray = qMinPair(:);

    % Get unique values using the MATLAB "unique" function 
    uniqueValues = unique(flatArray);

    % Convert to a row vector
    qMinVec = uniqueValues';

    % Initialize the E0 output vector
    E0MinVec = zeros(1,length(qMinVec)); 

    % Compare all values of qInputVec to qOutputVec. If they match,
    % index out the corresponding E0 value and store it in E0OutputVec. 
    for j=1:length(qMinVec)    
        for i=1:length(qVec)
            if qMinVec(j) == qVec(i)
                E0MinVec(j) = E0Vec(i);
            end
        end
    end

    %% %%%%%%%%% Creating and Displaying Graphs %%%%%%%%%%%%

    %% Plot 1: All Intersections

    % Initialize the figure
    figure
    grid on % add grid
    hold on
    
    % Define the Ef value range in order to plot the lines
    Ef_range = linspace(-1, 2, 1000);
    
    % Loop through all the charge states and plot each line
    for i = 1:length(qVec)
        E = E0Vec(i) + qVec(i)*Ef_range; % calculate the formation energy. 
        plot(Ef_range, E, 'DisplayName', sprintf('q = %d, E0 = %.2f eV', qVec(i), E0Vec(i))) % Plot the line and label
    end
    
    % Add Labels and Legend for the lines
    xlabel('Fermi Energy (E_f) [eV]');
    ylabel('Formation Energy (E) [eV]');
    title('Formation Energy vs. Fermi Energy for Different Charge States');
    legend('Location', 'best', 'Interpreter', 'none'); 
    
    % Plot the Intersection Points and Add Labels
    for j = 1:length(CTLVec)
        % Plot each intersection point without adding it to the legend
        plot(CTLVec(j), EVec(j), 'ro', 'MarkerSize', 8, 'HandleVisibility', 'off');
    
        % Add a label near the intersection point
        text(CTLVec(j), EVec(j), sprintf('([%0.1f,%0.1f], CTL = %0.2f)', ...
            qPair(j,1), qPair(j,2), CTLVec(j)), ...
            'VerticalAlignment', 'cap', 'HorizontalAlignment', 'left', ...
            'FontSize', 8, 'Color', 'black');
    end
    
    hold off

    %% Plot 2: Lowest Energy Lines
    
    % Initialize figure
    figure
    grid on % add grid
    hold on 
    
    % Loop through all the minimum charge states and plot each line
    for k = 1:length(qMinVec)
        E = E0MinVec(k) + qMinVec(k)*Ef_range; % calculate the formation energy. 
        plot(Ef_range, E, 'DisplayName', sprintf('q = %d, E0 = %.2f eV', qMinVec(k), E0MinVec(k))) % Plot the line and label
    end
    
    % Add Labels and Legend for the lines
    xlabel('Fermi Energy (E_f) [eV]');
    ylabel('Formation Energy (E) [eV]');
    title('Formation Energy vs. Fermi Energy for Minimum Charge States');
    legend('Location', 'best', 'Interpreter', 'none'); 
    
    % Plot the Intersection Points and Add Labels
    for p = 1:length(CTLMinVec)
        % Plot each intersection point without adding it to the legend
        plot(CTLMinVec(p), EMinVec(p), 'ro', 'MarkerSize', 8, 'HandleVisibility', 'off');
    
        % Add a label near the intersection point
        text(CTLMinVec(p), EMinVec(p), sprintf('([%0.1f,%0.1f], CTL = %0.2f)', ...
            qMinPair(p,1), qMinPair(p,2), CTLMinVec(p)), ...
            'VerticalAlignment', 'cap', 'HorizontalAlignment', 'left', ...
            'FontSize', 8, 'Color', 'black');
    end
    
    hold off

    %% Plot 3: Piecewise Lowest Energy Lines 
    % Initialize figure
    figure
    grid on % add grid
    hold on
    
    % Initialize variables to track the lowest formation energy and corresponding charge states
    lowestEnergy = inf(1, length(Ef_range)); % Initialize to very high values for comparison
    lowestChargeState = zeros(1, length(Ef_range)); % Store the charge state that gives the lowest energy
    
    % Loop through all charge states
    for i = 1:length(qVec)
        % Calculate the formation energy for the current charge state across all Ef values
        E = E0Vec(i) + qVec(i) * Ef_range; 
        
        % Compare with current lowest energy for each value of Ef
        for j = 1:length(Ef_range)
            if E(j) < lowestEnergy(j)
                lowestEnergy(j) = E(j); % Update to the new lowest energy
                lowestChargeState(j) = qVec(i); % Update to the corresponding charge state
            end
        end
    end
    
    % Plot the piecewise lowest energy line
    plot(Ef_range, lowestEnergy, 'blue', 'LineWidth', 1, 'DisplayName', 'Lowest Energy Line');
    
    % Add Labels and Legend for the lines
    xlabel('Fermi Energy (E_f) [eV]');
    ylabel('Formation Energy (E) [eV]');
    title('Piecewise Formation Energy vs. Fermi Energy for Minimum Charge States');
    legend('Location', 'northwest', 'Interpreter', 'none'); 
    
    % Highlight the transitions between different charge states
    % Loop through the range to identify changes in the lowest charge state
    p = length(CTLMinVec); % Initialize the counter variable for the text
    for j = 2:length(Ef_range)
        if lowestChargeState(j) ~= lowestChargeState(j-1)
            disp(p)
            % Plot a marker at the transition point
            plot(Ef_range(j), lowestEnergy(j), 'ro', 'MarkerSize', 8, 'HandleVisibility', 'off');
            % Add a label near the transition point to indicate the new charge state
            text(Ef_range(j), lowestEnergy(j), sprintf('([%0.1f,%0.1f], CTL = %0.2f)', ...
                qMinPair(p,1), qMinPair(p,2), CTLMinVec(p)), ...
            'VerticalAlignment', 'cap', 'HorizontalAlignment', 'left', ...
            'FontSize', 8, 'Color', 'black');
            p = p - 1; % update p
        end
    end
    
    hold off

    %% Transpose qMinVec and E0MinVec
    qMinVec = qMinVec(:); 
    E0MinVec = E0MinVec(:); 

end % end of the function
% function takes in a solution and defects variable, and for the specified
% element number, spits out the fractional contribution to that element's
% total concentration by chargestate and by defect.  This only makes sense
% for the impurities

function [defects_with_element_fraction, chargestates_with_element_fraction, defects_with_element_names, chargestates_with_element_names] = contains_element(solution,defects,element_num)
    
    chargestates_contains_index = find(defects.cs_dm(:,element_num)~=0);
    defects_contains_index = unique(defects.cs_ID(chargestates_contains_index));

    num_chargestates_contains = numel(chargestates_contains_index);
    num_defects_contains = numel(defects_contains_index);
    num_Ts = numel(solution.T_equilibrium);

    defects_with_element_names = defects.defect_names(defects_contains_index);
    chargestates_with_element_names = defects.chargestate_names(chargestates_contains_index);

    chargestates_with_element = solution.chargestates(:,chargestates_contains_index);
    defects_with_element = solution.defects(:,defects_contains_index);

    element_sum = sum((ones(num_Ts,1) * defects.cs_dm(chargestates_contains_index,element_num)')  .* chargestates_with_element , 2);   % this gives the scalar total amount of the element for each T
    chargestates_with_element_fraction = chargestates_with_element ./ (element_sum*ones(1,num_chargestates_contains));
    
    defects_with_element_fraction = zeros(num_Ts,num_defects_contains);


    for i = 1:num_defects_contains
        cs_lo = defects.cs_indices_lo(defects_contains_index(i));
        cs_hi = defects.cs_indices_hi(defects_contains_index(i));
        atoms_in_defect_i = sum((ones(num_Ts,1) * defects.cs_dm(cs_lo:cs_hi,element_num)') .* solution.chargestates(:,cs_lo:cs_hi),2);
        defects_with_element_fraction(:,i) = atoms_in_defect_i ./ element_sum;
    end

end
